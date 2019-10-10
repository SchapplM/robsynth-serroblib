% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:54
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
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:54
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
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(9);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = t47 * t53 - t48 * t54;
	t45 = t47 * t54 + t48 * t53;
	t44 = t47 * t52 + t48 * t55;
	t43 = t47 * t55 - t48 * t52;
	t1 = [t46, t43, 0, 0, 0, 0; -t44, -t45, 0, 0, 0, 0; 0, -qJD(2) * t47, 0, 0, 0, 0; t45, t44, 0, 0, 0, 0; t43, t46, 0, 0, 0, 0; 0, -qJD(2) * t48, 0, 0, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t180 = sin(qJ(1));
	t185 = qJD(1) * t180;
	t181 = cos(qJ(1));
	t184 = qJD(1) * t181;
	t183 = qJD(2) * t180;
	t182 = qJD(2) * t181;
	t179 = qJ(2) + pkin(9);
	t178 = cos(t179);
	t177 = sin(t179);
	t176 = -t177 * t183 + t178 * t184;
	t175 = t177 * t184 + t178 * t183;
	t174 = t177 * t182 + t178 * t185;
	t173 = -t177 * t185 + t178 * t182;
	t1 = [-t185, 0, 0, 0, 0, 0; t184, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t176, t173, 0, 0, 0, 0; t174, t175, 0, 0, 0, 0; 0, qJD(2) * t177, 0, 0, 0, 0; -t175, -t174, 0, 0, 0, 0; t173, t176, 0, 0, 0, 0; 0, qJD(2) * t178, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:55
	% EndTime: 2019-10-10 09:19:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->15), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
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
	t202 = qJ(2) + pkin(9);
	t201 = cos(t202);
	t213 = qJD(2) * t201;
	t212 = qJD(2) * t205;
	t211 = qJD(2) * t206;
	t210 = t201 * t212;
	t209 = t201 * t211;
	t200 = sin(t202);
	t208 = -t200 * t212 + t201 * t214;
	t207 = -t200 * t211 - t201 * t215;
	t1 = [-t203 * t210 + (-t200 * t218 - t217) * qJD(1), t207 * t203, 0, 0, 0, 0; t203 * t209 + (-t200 * t219 + t216) * qJD(1), t208 * t203, 0, 0, 0, 0; 0, t203 * t213, 0, 0, 0, 0; -t204 * t210 + (-t200 * t216 + t219) * qJD(1), t207 * t204, 0, 0, 0, 0; t204 * t209 + (-t200 * t217 - t218) * qJD(1), t208 * t204, 0, 0, 0, 0; 0, t204 * t213, 0, 0, 0, 0; -t208, t200 * t215 - t209, 0, 0, 0, 0; t207, -t200 * t214 - t210, 0, 0, 0, 0; 0, -qJD(2) * t200, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:55
	% EndTime: 2019-10-10 09:19:55
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (162->26), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->33)
	t279 = sin(qJ(1));
	t278 = qJ(2) + pkin(9);
	t274 = sin(t278);
	t292 = qJD(6) * t274;
	t286 = qJD(1) + t292;
	t300 = t279 * t286;
	t280 = cos(qJ(1));
	t299 = t280 * t286;
	t298 = qJD(1) * t279;
	t297 = qJD(1) * t280;
	t296 = qJD(2) * t274;
	t276 = cos(t278);
	t295 = qJD(2) * t276;
	t294 = qJD(2) * t279;
	t293 = qJD(2) * t280;
	t291 = qJD(6) * t276;
	t277 = pkin(10) + qJ(6);
	t273 = sin(t277);
	t290 = t273 * t291;
	t275 = cos(t277);
	t289 = t275 * t291;
	t288 = t276 * t294;
	t287 = t276 * t293;
	t285 = -qJD(1) * t274 - qJD(6);
	t284 = -t274 * t294 + t276 * t297;
	t283 = -t274 * t293 - t276 * t298;
	t282 = t285 * t280 - t288;
	t281 = t285 * t279 + t287;
	t272 = t281 * t273 + t275 * t299;
	t271 = -t273 * t299 + t281 * t275;
	t270 = t282 * t273 - t275 * t300;
	t269 = t273 * t300 + t282 * t275;
	t1 = [t270, t283 * t273 + t280 * t289, 0, 0, 0, t271; t272, t284 * t273 + t279 * t289, 0, 0, 0, -t269; 0, t273 * t295 + t275 * t292, 0, 0, 0, t275 * t296 + t290; t269, t283 * t275 - t280 * t290, 0, 0, 0, -t272; t271, t284 * t275 - t279 * t290, 0, 0, 0, t270; 0, -t273 * t292 + t275 * t295, 0, 0, 0, -t273 * t296 + t289; -t284, t274 * t298 - t287, 0, 0, 0, 0; t283, -t274 * t297 - t288, 0, 0, 0, 0; 0, -t296, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end