% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t24 = qJD(1) * sin(qJ(1));
	t23 = qJD(1) * cos(qJ(1));
	t20 = cos(pkin(9));
	t19 = sin(pkin(9));
	t1 = [-t19 * t24, 0, 0, 0, 0, 0; t19 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t20 * t24, 0, 0, 0, 0, 0; t20 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t23, 0, 0, 0, 0, 0; -t24, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t44 = sin(qJ(1));
	t49 = qJD(1) * t44;
	t45 = cos(qJ(1));
	t48 = qJD(1) * t45;
	t47 = qJD(4) * t44;
	t46 = qJD(4) * t45;
	t43 = pkin(9) + qJ(4);
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = -t41 * t47 + t42 * t48;
	t39 = t41 * t48 + t42 * t47;
	t38 = t41 * t46 + t42 * t49;
	t37 = -t41 * t49 + t42 * t46;
	t1 = [t37, 0, 0, t40, 0, 0; t39, 0, 0, t38, 0, 0; 0, 0, 0, -qJD(4) * t42, 0, 0; -t38, 0, 0, -t39, 0, 0; t40, 0, 0, t37, 0, 0; 0, 0, 0, qJD(4) * t41, 0, 0; -t48, 0, 0, 0, 0, 0; -t49, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:25
	% EndTime: 2019-10-09 23:44:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->18), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t196 = sin(pkin(10));
	t198 = sin(qJ(1));
	t212 = t196 * t198;
	t199 = cos(qJ(1));
	t211 = t196 * t199;
	t197 = cos(pkin(10));
	t210 = t197 * t198;
	t209 = t197 * t199;
	t208 = qJD(1) * t198;
	t207 = qJD(1) * t199;
	t195 = pkin(9) + qJ(4);
	t194 = cos(t195);
	t206 = qJD(4) * t194;
	t205 = qJD(4) * t198;
	t204 = qJD(4) * t199;
	t203 = t194 * t205;
	t202 = t194 * t204;
	t193 = sin(t195);
	t201 = -t193 * t205 + t194 * t207;
	t200 = t193 * t204 + t194 * t208;
	t1 = [t197 * t202 + (-t193 * t210 - t211) * qJD(1), 0, 0, t201 * t197, 0, 0; t197 * t203 + (t193 * t209 - t212) * qJD(1), 0, 0, t200 * t197, 0, 0; 0, 0, 0, -t197 * t206, 0, 0; -t196 * t202 + (t193 * t212 - t209) * qJD(1), 0, 0, -t201 * t196, 0, 0; -t196 * t203 + (-t193 * t211 - t210) * qJD(1), 0, 0, -t200 * t196, 0, 0; 0, 0, 0, t196 * t206, 0, 0; t200, 0, 0, t193 * t207 + t203, 0, 0; -t201, 0, 0, t193 * t208 - t202, 0, 0; 0, 0, 0, -qJD(4) * t193, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:26
	% EndTime: 2019-10-09 23:44:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (162->26), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->33)
	t268 = sin(qJ(1));
	t267 = pkin(9) + qJ(4);
	t263 = sin(t267);
	t274 = qJD(1) * t263 + qJD(6);
	t265 = cos(t267);
	t269 = cos(qJ(1));
	t282 = qJD(4) * t269;
	t276 = t265 * t282;
	t289 = t274 * t268 - t276;
	t283 = qJD(4) * t268;
	t277 = t265 * t283;
	t288 = t274 * t269 + t277;
	t287 = qJD(1) * t268;
	t286 = qJD(1) * t269;
	t285 = qJD(4) * t263;
	t284 = qJD(4) * t265;
	t281 = qJD(6) * t263;
	t280 = qJD(6) * t265;
	t266 = pkin(10) + qJ(6);
	t262 = sin(t266);
	t279 = t262 * t280;
	t264 = cos(t266);
	t278 = t264 * t280;
	t275 = -qJD(1) - t281;
	t273 = t275 * t268;
	t272 = t275 * t269;
	t271 = -t263 * t283 + t265 * t286;
	t270 = t263 * t282 + t265 * t287;
	t261 = t262 * t273 + t288 * t264;
	t260 = -t288 * t262 + t264 * t273;
	t259 = t262 * t272 - t289 * t264;
	t258 = t289 * t262 + t264 * t272;
	t1 = [t259, 0, 0, t271 * t264 - t268 * t279, 0, t260; t261, 0, 0, t270 * t264 + t269 * t279, 0, -t258; 0, 0, 0, t262 * t281 - t264 * t284, 0, t262 * t285 - t278; t258, 0, 0, -t271 * t262 - t268 * t278, 0, -t261; t260, 0, 0, -t270 * t262 + t269 * t278, 0, t259; 0, 0, 0, t262 * t284 + t264 * t281, 0, t264 * t285 + t279; t270, 0, 0, t263 * t286 + t277, 0, 0; -t271, 0, 0, t263 * t287 - t276, 0, 0; 0, 0, 0, -t285, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end