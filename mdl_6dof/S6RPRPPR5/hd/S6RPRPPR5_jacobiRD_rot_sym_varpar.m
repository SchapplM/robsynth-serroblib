% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
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
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(9));
	t17 = sin(pkin(9));
	t1 = [-t18 * t21, 0, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:18
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(9) + qJ(3);
	t41 = cos(t42);
	t40 = sin(t42);
	t39 = t40 * t46 - t41 * t47;
	t38 = t40 * t47 + t41 * t46;
	t37 = t40 * t45 + t41 * t48;
	t36 = t40 * t48 - t41 * t45;
	t1 = [t39, 0, t36, 0, 0, 0; -t37, 0, -t38, 0, 0, 0; 0, 0, -qJD(3) * t40, 0, 0, 0; t38, 0, t37, 0, 0, 0; t36, 0, t39, 0, 0, 0; 0, 0, -qJD(3) * t41, 0, 0, 0; -t48, 0, 0, 0, 0, 0; t47, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:18
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t170 = sin(qJ(1));
	t175 = qJD(1) * t170;
	t171 = cos(qJ(1));
	t174 = qJD(1) * t171;
	t173 = qJD(3) * t170;
	t172 = qJD(3) * t171;
	t169 = pkin(9) + qJ(3);
	t168 = cos(t169);
	t167 = sin(t169);
	t166 = -t167 * t173 + t168 * t174;
	t165 = t167 * t174 + t168 * t173;
	t164 = t167 * t172 + t168 * t175;
	t163 = -t167 * t175 + t168 * t172;
	t1 = [-t175, 0, 0, 0, 0, 0; t174, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t166, 0, t163, 0, 0, 0; t164, 0, t165, 0, 0, 0; 0, 0, qJD(3) * t167, 0, 0, 0; -t165, 0, -t164, 0, 0, 0; t163, 0, t166, 0, 0, 0; 0, 0, qJD(3) * t168, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:19
	% EndTime: 2019-10-10 00:22:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->15), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t195 = sin(pkin(10));
	t197 = sin(qJ(1));
	t211 = t195 * t197;
	t198 = cos(qJ(1));
	t210 = t195 * t198;
	t196 = cos(pkin(10));
	t209 = t196 * t197;
	t208 = t196 * t198;
	t207 = qJD(1) * t197;
	t206 = qJD(1) * t198;
	t194 = pkin(9) + qJ(3);
	t193 = cos(t194);
	t205 = qJD(3) * t193;
	t204 = qJD(3) * t197;
	t203 = qJD(3) * t198;
	t202 = t193 * t204;
	t201 = t193 * t203;
	t192 = sin(t194);
	t200 = -t192 * t204 + t193 * t206;
	t199 = -t192 * t203 - t193 * t207;
	t1 = [-t195 * t202 + (-t192 * t210 - t209) * qJD(1), 0, t199 * t195, 0, 0, 0; t195 * t201 + (-t192 * t211 + t208) * qJD(1), 0, t200 * t195, 0, 0, 0; 0, 0, t195 * t205, 0, 0, 0; -t196 * t202 + (-t192 * t208 + t211) * qJD(1), 0, t199 * t196, 0, 0, 0; t196 * t201 + (-t192 * t209 - t210) * qJD(1), 0, t200 * t196, 0, 0, 0; 0, 0, t196 * t205, 0, 0, 0; -t200, 0, t192 * t207 - t201, 0, 0, 0; t199, 0, -t192 * t206 - t202, 0, 0, 0; 0, 0, -qJD(3) * t192, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:19
	% EndTime: 2019-10-10 00:22:19
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (162->26), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->33)
	t272 = sin(qJ(1));
	t271 = pkin(9) + qJ(3);
	t267 = sin(t271);
	t285 = qJD(6) * t267;
	t279 = qJD(1) + t285;
	t293 = t272 * t279;
	t273 = cos(qJ(1));
	t292 = t273 * t279;
	t291 = qJD(1) * t272;
	t290 = qJD(1) * t273;
	t289 = qJD(3) * t267;
	t269 = cos(t271);
	t288 = qJD(3) * t269;
	t287 = qJD(3) * t272;
	t286 = qJD(3) * t273;
	t284 = qJD(6) * t269;
	t270 = pkin(10) + qJ(6);
	t266 = sin(t270);
	t283 = t266 * t284;
	t268 = cos(t270);
	t282 = t268 * t284;
	t281 = t269 * t287;
	t280 = t269 * t286;
	t278 = -qJD(1) * t267 - qJD(6);
	t277 = -t267 * t287 + t269 * t290;
	t276 = -t267 * t286 - t269 * t291;
	t275 = t278 * t273 - t281;
	t274 = t278 * t272 + t280;
	t265 = t274 * t266 + t268 * t292;
	t264 = -t266 * t292 + t274 * t268;
	t263 = t275 * t266 - t268 * t293;
	t262 = t266 * t293 + t275 * t268;
	t1 = [t263, 0, t276 * t266 + t273 * t282, 0, 0, t264; t265, 0, t277 * t266 + t272 * t282, 0, 0, -t262; 0, 0, t266 * t288 + t268 * t285, 0, 0, t268 * t289 + t283; t262, 0, t276 * t268 - t273 * t283, 0, 0, -t265; t264, 0, t277 * t268 - t272 * t283, 0, 0, t263; 0, 0, -t266 * t285 + t268 * t288, 0, 0, -t266 * t289 + t282; -t277, 0, t267 * t291 - t280, 0, 0, 0; t276, 0, -t267 * t290 - t281, 0, 0, 0; 0, 0, -t289, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end