% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
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
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t16 = qJD(1) * sin(qJ(1));
	t15 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t15, 0, 0, 0, 0, 0; -t16, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:05
	% EndTime: 2019-10-09 23:41:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(4));
	t40 = qJD(4) * t34;
	t36 = cos(qJ(4));
	t39 = qJD(4) * t36;
	t38 = qJD(4) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = -t34 * t41 - t35 * t39;
	t31 = -t34 * t38 - t36 * t42;
	t30 = t34 * t42 - t36 * t38;
	t1 = [t32, 0, 0, t31, 0, 0; -t30, 0, 0, t33, 0, 0; 0, 0, 0, -t39, 0, 0; -t33, 0, 0, t30, 0, 0; t31, 0, 0, t32, 0, 0; 0, 0, 0, t40, 0, 0; t42, 0, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:06
	% EndTime: 2019-10-09 23:41:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (18->16), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t185 = sin(qJ(4));
	t186 = sin(qJ(1));
	t199 = t185 * t186;
	t188 = cos(qJ(1));
	t198 = t185 * t188;
	t197 = qJD(1) * t186;
	t196 = qJD(1) * t188;
	t195 = qJD(4) * t185;
	t187 = cos(qJ(4));
	t194 = qJD(4) * t187;
	t193 = qJD(4) * t188;
	t192 = t186 * t194;
	t191 = t187 * t193;
	t190 = -t186 * t195 + t187 * t196;
	t189 = t185 * t193 + t187 * t197;
	t184 = cos(pkin(9));
	t183 = sin(pkin(9));
	t1 = [-t184 * t192 + (t183 * t186 - t184 * t198) * qJD(1), 0, 0, -t189 * t184, 0, 0; t184 * t191 + (-t183 * t188 - t184 * t199) * qJD(1), 0, 0, t190 * t184, 0, 0; 0, 0, 0, -t184 * t194, 0, 0; t183 * t192 + (t183 * t198 + t184 * t186) * qJD(1), 0, 0, t189 * t183, 0, 0; -t183 * t191 + (t183 * t199 - t184 * t188) * qJD(1), 0, 0, -t190 * t183, 0, 0; 0, 0, 0, t183 * t194, 0, 0; t190, 0, 0, -t185 * t197 + t191, 0, 0; t189, 0, 0, t185 * t196 + t192, 0, 0; 0, 0, 0, -t195, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:06
	% EndTime: 2019-10-09 23:41:06
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (109->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t255 = cos(qJ(1));
	t252 = sin(qJ(4));
	t267 = qJD(6) * t252;
	t261 = qJD(1) + t267;
	t274 = t255 * t261;
	t260 = qJD(1) * t252 + qJD(6);
	t253 = sin(qJ(1));
	t254 = cos(qJ(4));
	t269 = qJD(4) * t254;
	t265 = t253 * t269;
	t273 = t260 * t255 + t265;
	t272 = qJD(1) * t253;
	t271 = qJD(1) * t255;
	t270 = qJD(4) * t252;
	t268 = qJD(4) * t255;
	t266 = qJD(6) * t254;
	t251 = pkin(9) + qJ(6);
	t249 = sin(t251);
	t264 = t249 * t266;
	t250 = cos(t251);
	t263 = t250 * t266;
	t262 = t254 * t268;
	t259 = t261 * t253;
	t258 = -t253 * t270 + t254 * t271;
	t257 = t252 * t268 + t254 * t272;
	t256 = t260 * t253 - t262;
	t248 = t249 * t259 - t273 * t250;
	t247 = t273 * t249 + t250 * t259;
	t246 = t249 * t274 + t256 * t250;
	t245 = t256 * t249 - t250 * t274;
	t1 = [t248, 0, 0, -t257 * t250 - t255 * t264, 0, t245; -t246, 0, 0, t258 * t250 - t253 * t264, 0, -t247; 0, 0, 0, t249 * t267 - t250 * t269, 0, t249 * t270 - t263; t247, 0, 0, t257 * t249 - t255 * t263, 0, t246; t245, 0, 0, -t258 * t249 - t253 * t263, 0, t248; 0, 0, 0, t249 * t269 + t250 * t267, 0, t250 * t270 + t264; t258, 0, 0, -t252 * t272 + t262, 0, 0; t257, 0, 0, t252 * t271 + t265, 0, 0; 0, 0, 0, -t270, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end