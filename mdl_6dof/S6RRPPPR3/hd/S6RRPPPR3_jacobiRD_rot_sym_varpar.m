% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:21
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:37
	% EndTime: 2019-10-10 09:21:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t156 = sin(qJ(1));
	t163 = qJD(1) * t156;
	t158 = cos(qJ(1));
	t162 = qJD(1) * t158;
	t155 = sin(qJ(2));
	t161 = qJD(2) * t155;
	t157 = cos(qJ(2));
	t160 = qJD(2) * t157;
	t159 = qJD(2) * t158;
	t154 = -t156 * t161 + t157 * t162;
	t153 = -t155 * t162 - t156 * t160;
	t152 = -t155 * t159 - t157 * t163;
	t151 = t155 * t163 - t157 * t159;
	t1 = [-t154, t151, 0, 0, 0, 0; t152, t153, 0, 0, 0, 0; 0, -t161, 0, 0, 0, 0; -t163, 0, 0, 0, 0, 0; t162, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t153, t152, 0, 0, 0, 0; -t151, t154, 0, 0, 0, 0; 0, t160, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:37
	% EndTime: 2019-10-10 09:21:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(1));
	t48 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t47 = qJD(1) * t43;
	t40 = sin(qJ(2));
	t46 = qJD(2) * t40;
	t42 = cos(qJ(2));
	t45 = qJD(2) * t42;
	t44 = qJD(2) * t43;
	t39 = -t41 * t46 + t42 * t47;
	t38 = t40 * t47 + t41 * t45;
	t37 = t40 * t44 + t42 * t48;
	t36 = -t40 * t48 + t42 * t44;
	t1 = [-t38, -t37, 0, 0, 0, 0; t36, t39, 0, 0, 0, 0; 0, t45, 0, 0, 0, 0; t39, t36, 0, 0, 0, 0; t37, t38, 0, 0, 0, 0; 0, t46, 0, 0, 0, 0; t48, 0, 0, 0, 0, 0; -t47, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:38
	% EndTime: 2019-10-10 09:21:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (18->18), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t187 = sin(qJ(2));
	t188 = sin(qJ(1));
	t201 = t187 * t188;
	t190 = cos(qJ(1));
	t200 = t187 * t190;
	t199 = qJD(1) * t188;
	t198 = qJD(1) * t190;
	t197 = qJD(2) * t187;
	t189 = cos(qJ(2));
	t196 = qJD(2) * t189;
	t195 = qJD(2) * t190;
	t194 = t188 * t196;
	t193 = t189 * t195;
	t192 = -t188 * t197 + t189 * t198;
	t191 = t187 * t195 + t189 * t199;
	t186 = cos(pkin(9));
	t185 = sin(pkin(9));
	t1 = [-t186 * t194 + (t185 * t188 - t186 * t200) * qJD(1), -t191 * t186, 0, 0, 0, 0; t186 * t193 + (-t185 * t190 - t186 * t201) * qJD(1), t192 * t186, 0, 0, 0, 0; 0, t186 * t196, 0, 0, 0, 0; t185 * t194 + (t185 * t200 + t186 * t188) * qJD(1), t191 * t185, 0, 0, 0, 0; -t185 * t193 + (t185 * t201 - t186 * t190) * qJD(1), -t192 * t185, 0, 0, 0, 0; 0, -t185 * t196, 0, 0, 0, 0; -t192, t187 * t199 - t193, 0, 0, 0, 0; -t191, -t187 * t198 - t194, 0, 0, 0, 0; 0, -t197, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:38
	% EndTime: 2019-10-10 09:21:38
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (109->26), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t257 = cos(qJ(1));
	t254 = sin(qJ(2));
	t269 = qJD(6) * t254;
	t263 = qJD(1) + t269;
	t276 = t257 * t263;
	t262 = qJD(1) * t254 + qJD(6);
	t255 = sin(qJ(1));
	t256 = cos(qJ(2));
	t271 = qJD(2) * t256;
	t267 = t255 * t271;
	t275 = t262 * t257 + t267;
	t274 = qJD(1) * t255;
	t273 = qJD(1) * t257;
	t272 = qJD(2) * t254;
	t270 = qJD(2) * t257;
	t268 = qJD(6) * t256;
	t253 = pkin(9) + qJ(6);
	t251 = sin(t253);
	t266 = t251 * t268;
	t252 = cos(t253);
	t265 = t252 * t268;
	t264 = t256 * t270;
	t261 = t263 * t255;
	t260 = -t255 * t272 + t256 * t273;
	t259 = t254 * t270 + t256 * t274;
	t258 = t262 * t255 - t264;
	t250 = t251 * t261 - t275 * t252;
	t249 = t275 * t251 + t252 * t261;
	t248 = t251 * t276 + t258 * t252;
	t247 = t258 * t251 - t252 * t276;
	t1 = [t250, -t259 * t252 - t257 * t266, 0, 0, 0, t247; -t248, t260 * t252 - t255 * t266, 0, 0, 0, -t249; 0, -t251 * t269 + t252 * t271, 0, 0, 0, -t251 * t272 + t265; t249, t259 * t251 - t257 * t265, 0, 0, 0, t248; t247, -t260 * t251 - t255 * t265, 0, 0, 0, t250; 0, -t251 * t271 - t252 * t269, 0, 0, 0, -t252 * t272 - t266; -t260, t254 * t274 - t264, 0, 0, 0, 0; -t259, -t254 * t273 - t267, 0, 0, 0, 0; 0, -t272, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end