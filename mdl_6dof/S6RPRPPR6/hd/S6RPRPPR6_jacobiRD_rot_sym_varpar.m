% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:00
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
	% StartTime: 2019-10-10 00:24:00
	% EndTime: 2019-10-10 00:24:01
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
	% StartTime: 2019-10-10 00:24:01
	% EndTime: 2019-10-10 00:24:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0, 0, 0; t32, 0, t31, 0, 0, 0; 0, 0, -t39, 0, 0, 0; -t31, 0, -t32, 0, 0, 0; t33, 0, t30, 0, 0, 0; 0, 0, t40, 0, 0, 0; -t41, 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:01
	% EndTime: 2019-10-10 00:24:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t47 = sin(qJ(1));
	t52 = qJD(1) * t47;
	t48 = cos(qJ(1));
	t51 = qJD(1) * t48;
	t50 = qJD(3) * t47;
	t49 = qJD(3) * t48;
	t46 = qJ(3) + pkin(9);
	t45 = cos(t46);
	t44 = sin(t46);
	t43 = -t44 * t50 + t45 * t51;
	t42 = t44 * t51 + t45 * t50;
	t41 = t44 * t49 + t45 * t52;
	t40 = -t44 * t52 + t45 * t49;
	t1 = [t40, 0, t43, 0, 0, 0; t42, 0, t41, 0, 0, 0; 0, 0, -qJD(3) * t45, 0, 0, 0; -t41, 0, -t42, 0, 0, 0; t43, 0, t40, 0, 0, 0; 0, 0, qJD(3) * t44, 0, 0, 0; -t51, 0, 0, 0, 0, 0; -t52, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:02
	% EndTime: 2019-10-10 00:24:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->18), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t203 = sin(pkin(10));
	t205 = sin(qJ(1));
	t219 = t203 * t205;
	t206 = cos(qJ(1));
	t218 = t203 * t206;
	t204 = cos(pkin(10));
	t217 = t204 * t205;
	t216 = t204 * t206;
	t215 = qJD(1) * t205;
	t214 = qJD(1) * t206;
	t202 = qJ(3) + pkin(9);
	t201 = cos(t202);
	t213 = qJD(3) * t201;
	t212 = qJD(3) * t205;
	t211 = qJD(3) * t206;
	t210 = t201 * t212;
	t209 = t201 * t211;
	t200 = sin(t202);
	t208 = -t200 * t212 + t201 * t214;
	t207 = t200 * t211 + t201 * t215;
	t1 = [t204 * t209 + (-t200 * t217 - t218) * qJD(1), 0, t208 * t204, 0, 0, 0; t204 * t210 + (t200 * t216 - t219) * qJD(1), 0, t207 * t204, 0, 0, 0; 0, 0, -t204 * t213, 0, 0, 0; -t203 * t209 + (t200 * t219 - t216) * qJD(1), 0, -t208 * t203, 0, 0, 0; -t203 * t210 + (-t200 * t218 - t217) * qJD(1), 0, -t207 * t203, 0, 0, 0; 0, 0, t203 * t213, 0, 0, 0; t207, 0, t200 * t214 + t210, 0, 0, 0; -t208, 0, t200 * t215 - t209, 0, 0, 0; 0, 0, -qJD(3) * t200, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:24:02
	% EndTime: 2019-10-10 00:24:02
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (162->26), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->33)
	t274 = sin(qJ(1));
	t273 = qJ(3) + pkin(9);
	t269 = sin(t273);
	t280 = qJD(1) * t269 + qJD(6);
	t271 = cos(t273);
	t275 = cos(qJ(1));
	t288 = qJD(3) * t275;
	t282 = t271 * t288;
	t295 = t280 * t274 - t282;
	t289 = qJD(3) * t274;
	t283 = t271 * t289;
	t294 = t280 * t275 + t283;
	t293 = qJD(1) * t274;
	t292 = qJD(1) * t275;
	t291 = qJD(3) * t269;
	t290 = qJD(3) * t271;
	t287 = qJD(6) * t269;
	t286 = qJD(6) * t271;
	t272 = pkin(10) + qJ(6);
	t268 = sin(t272);
	t285 = t268 * t286;
	t270 = cos(t272);
	t284 = t270 * t286;
	t281 = -qJD(1) - t287;
	t279 = t281 * t274;
	t278 = t281 * t275;
	t277 = -t269 * t289 + t271 * t292;
	t276 = t269 * t288 + t271 * t293;
	t267 = t268 * t279 + t294 * t270;
	t266 = -t294 * t268 + t270 * t279;
	t265 = t268 * t278 - t295 * t270;
	t264 = t295 * t268 + t270 * t278;
	t1 = [t265, 0, t277 * t270 - t274 * t285, 0, 0, t266; t267, 0, t276 * t270 + t275 * t285, 0, 0, -t264; 0, 0, t268 * t287 - t270 * t290, 0, 0, t268 * t291 - t284; t264, 0, -t277 * t268 - t274 * t284, 0, 0, -t267; t266, 0, -t276 * t268 + t275 * t284, 0, 0, t265; 0, 0, t268 * t290 + t270 * t287, 0, 0, t270 * t291 + t285; t276, 0, t269 * t292 + t283, 0, 0, 0; -t277, 0, t269 * t293 - t282, 0, 0, 0; 0, 0, -t291, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end