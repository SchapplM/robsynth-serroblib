% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
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
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(9);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->4), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(1) * sin(pkin(10));
	t25 = qJD(1) * cos(pkin(10));
	t22 = qJ(1) + pkin(9);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [-t21 * t25, 0, 0, 0, 0, 0; -t20 * t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t21 * t26, 0, 0, 0, 0, 0; t20 * t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -qJD(1) * t20, 0, 0, 0, 0, 0; qJD(1) * t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t50 = qJ(1) + pkin(9);
	t46 = sin(t50);
	t55 = qJD(1) * t46;
	t48 = cos(t50);
	t54 = qJD(1) * t48;
	t49 = pkin(10) + qJ(4);
	t45 = sin(t49);
	t53 = qJD(4) * t45;
	t47 = cos(t49);
	t52 = qJD(4) * t47;
	t51 = qJD(4) * t48;
	t44 = t46 * t53 - t47 * t54;
	t43 = t45 * t54 + t46 * t52;
	t42 = t45 * t51 + t47 * t55;
	t41 = t45 * t55 - t47 * t51;
	t1 = [t44, 0, 0, t41, 0, 0; -t42, 0, 0, -t43, 0, 0; 0, 0, 0, -t53, 0, 0; t43, 0, 0, t42, 0, 0; t41, 0, 0, t44, 0, 0; 0, 0, 0, -t52, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:01
	% EndTime: 2019-10-09 23:36:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t181 = qJ(1) + pkin(9);
	t177 = sin(t181);
	t186 = qJD(1) * t177;
	t179 = cos(t181);
	t185 = qJD(1) * t179;
	t180 = pkin(10) + qJ(4);
	t176 = sin(t180);
	t184 = qJD(4) * t176;
	t178 = cos(t180);
	t183 = qJD(4) * t178;
	t182 = qJD(4) * t179;
	t175 = -t177 * t184 + t178 * t185;
	t174 = t176 * t185 + t177 * t183;
	t173 = t176 * t182 + t178 * t186;
	t172 = -t176 * t186 + t178 * t182;
	t1 = [-t186, 0, 0, 0, 0, 0; t185, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t175, 0, 0, t172, 0, 0; t173, 0, 0, t174, 0, 0; 0, 0, 0, t184, 0, 0; -t174, 0, 0, -t173, 0, 0; t172, 0, 0, t175, 0, 0; 0, 0, 0, t183, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:02
	% EndTime: 2019-10-09 23:36:02
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (162->27), mult. (173->51), div. (0->0), fcn. (173->6), ass. (0->33)
	t258 = sin(qJ(6));
	t256 = pkin(10) + qJ(4);
	t252 = sin(t256);
	t265 = qJD(6) * t252 + qJD(1);
	t254 = cos(t256);
	t259 = cos(qJ(6));
	t272 = qJD(4) * t259;
	t267 = t254 * t272;
	t279 = t265 * t258 - t267;
	t273 = qJD(4) * t258;
	t268 = t254 * t273;
	t278 = t265 * t259 + t268;
	t257 = qJ(1) + pkin(9);
	t253 = sin(t257);
	t277 = qJD(1) * t253;
	t255 = cos(t257);
	t276 = qJD(1) * t255;
	t275 = qJD(4) * t253;
	t274 = qJD(4) * t255;
	t271 = qJD(6) * t258;
	t270 = qJD(6) * t259;
	t269 = t254 * t276;
	t266 = t252 * t274;
	t264 = -qJD(1) * t252 - qJD(6);
	t263 = t264 * t258;
	t262 = t264 * t259;
	t261 = t252 * t272 + t254 * t271;
	t260 = -t252 * t273 + t254 * t270;
	t251 = t253 * t263 + t278 * t255;
	t250 = t253 * t262 - t279 * t255;
	t249 = -t278 * t253 + t255 * t263;
	t248 = t279 * t253 + t255 * t262;
	t1 = [t249, 0, 0, -t258 * t266 + (t255 * t270 - t258 * t277) * t254, 0, t250; t251, 0, 0, t260 * t253 + t258 * t269, 0, -t248; 0, 0, 0, t252 * t270 + t268, 0, t261; t248, 0, 0, -t259 * t266 + (-t255 * t271 - t259 * t277) * t254, 0, -t251; t250, 0, 0, -t261 * t253 + t259 * t269, 0, t249; 0, 0, 0, -t252 * t271 + t267, 0, t260; t252 * t275 - t269, 0, 0, t252 * t277 - t254 * t274, 0, 0; -t254 * t277 - t266, 0, 0, -t252 * t276 - t254 * t275, 0, 0; 0, 0, 0, -qJD(4) * t252, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end