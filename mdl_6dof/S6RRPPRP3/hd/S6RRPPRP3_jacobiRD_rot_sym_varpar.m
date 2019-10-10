% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:30
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
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
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
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
	% StartTime: 2019-10-10 09:30:30
	% EndTime: 2019-10-10 09:30:30
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
	% StartTime: 2019-10-10 09:30:30
	% EndTime: 2019-10-10 09:30:30
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
	% StartTime: 2019-10-10 09:30:31
	% EndTime: 2019-10-10 09:30:31
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (49->25), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->32)
	t241 = cos(qJ(5));
	t243 = cos(qJ(1));
	t264 = t241 * t243;
	t240 = sin(qJ(1));
	t263 = qJD(1) * t240;
	t262 = qJD(1) * t243;
	t239 = sin(qJ(2));
	t261 = qJD(2) * t239;
	t242 = cos(qJ(2));
	t260 = qJD(2) * t242;
	t259 = qJD(2) * t243;
	t258 = qJD(5) * t239;
	t257 = qJD(5) * t242;
	t256 = t242 * t262;
	t255 = t241 * t260;
	t254 = t240 * t260;
	t238 = sin(qJ(5));
	t253 = t238 * t257;
	t252 = t241 * t257;
	t251 = t242 * t259;
	t250 = qJD(1) + t258;
	t249 = qJD(1) * t239 + qJD(5);
	t248 = t250 * t238;
	t247 = t240 * t261 - t256;
	t246 = t239 * t259 + t242 * t263;
	t245 = -t241 * t261 - t253;
	t244 = t249 * t240 - t251;
	t237 = -t249 * t264 + (t248 - t255) * t240;
	t236 = t250 * t241 * t240 + (t249 * t243 + t254) * t238;
	t235 = t244 * t241 + t243 * t248;
	t234 = t244 * t238 - t250 * t264;
	t1 = [t237, -t246 * t241 - t243 * t253, 0, 0, t234, 0; -t235, t245 * t240 + t241 * t256, 0, 0, -t236, 0; 0, -t238 * t258 + t255, 0, 0, -t238 * t261 + t252, 0; t236, t246 * t238 - t243 * t252, 0, 0, t235, 0; t234, t247 * t238 - t240 * t252, 0, 0, t237, 0; 0, -t238 * t260 - t241 * t258, 0, 0, t245, 0; t247, t239 * t263 - t251, 0, 0, 0, 0; -t246, -t239 * t262 - t254, 0, 0, 0, 0; 0, -t261, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:31
	% EndTime: 2019-10-10 09:30:31
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (49->25), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->32)
	t257 = cos(qJ(5));
	t259 = cos(qJ(1));
	t280 = t257 * t259;
	t256 = sin(qJ(1));
	t279 = qJD(1) * t256;
	t278 = qJD(1) * t259;
	t255 = sin(qJ(2));
	t277 = qJD(2) * t255;
	t258 = cos(qJ(2));
	t276 = qJD(2) * t258;
	t275 = qJD(2) * t259;
	t274 = qJD(5) * t255;
	t273 = qJD(5) * t258;
	t272 = t258 * t278;
	t271 = t257 * t276;
	t270 = t256 * t276;
	t254 = sin(qJ(5));
	t269 = t254 * t273;
	t268 = t257 * t273;
	t267 = t258 * t275;
	t266 = qJD(1) + t274;
	t265 = qJD(1) * t255 + qJD(5);
	t264 = t266 * t254;
	t263 = t256 * t277 - t272;
	t262 = t255 * t275 + t258 * t279;
	t261 = -t257 * t277 - t269;
	t260 = t265 * t256 - t267;
	t253 = -t265 * t280 + (t264 - t271) * t256;
	t252 = t266 * t257 * t256 + (t265 * t259 + t270) * t254;
	t251 = t260 * t257 + t259 * t264;
	t250 = t260 * t254 - t266 * t280;
	t1 = [t253, -t262 * t257 - t259 * t269, 0, 0, t250, 0; -t251, t261 * t256 + t257 * t272, 0, 0, -t252, 0; 0, -t254 * t274 + t271, 0, 0, -t254 * t277 + t268, 0; t252, t262 * t254 - t259 * t268, 0, 0, t251, 0; t250, t263 * t254 - t256 * t268, 0, 0, t253, 0; 0, -t254 * t276 - t257 * t274, 0, 0, t261, 0; t263, t255 * t279 - t267, 0, 0, 0, 0; -t262, -t255 * t278 - t270, 0, 0, 0, 0; 0, -t277, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end