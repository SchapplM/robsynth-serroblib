% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
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
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:09
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
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:09
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
	% StartTime: 2019-10-10 09:18:10
	% EndTime: 2019-10-10 09:18:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->18), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t204 = sin(pkin(10));
	t206 = sin(qJ(1));
	t220 = t204 * t206;
	t207 = cos(qJ(1));
	t219 = t204 * t207;
	t205 = cos(pkin(10));
	t218 = t205 * t206;
	t217 = t205 * t207;
	t216 = qJD(1) * t206;
	t215 = qJD(1) * t207;
	t203 = qJ(2) + pkin(9);
	t201 = sin(t203);
	t214 = qJD(2) * t201;
	t213 = qJD(2) * t206;
	t212 = qJD(2) * t207;
	t211 = t201 * t213;
	t210 = t201 * t212;
	t202 = cos(t203);
	t209 = t201 * t215 + t202 * t213;
	t208 = t201 * t216 - t202 * t212;
	t1 = [t205 * t211 + (-t202 * t217 - t220) * qJD(1), t208 * t205, 0, 0, 0, 0; -t205 * t210 + (-t202 * t218 + t219) * qJD(1), -t209 * t205, 0, 0, 0, 0; 0, -t205 * t214, 0, 0, 0, 0; -t204 * t211 + (t202 * t219 - t218) * qJD(1), -t208 * t204, 0, 0, 0, 0; t204 * t210 + (t202 * t220 + t217) * qJD(1), t209 * t204, 0, 0, 0, 0; 0, t204 * t214, 0, 0, 0, 0; -t209, -t202 * t216 - t210, 0, 0, 0, 0; -t208, t202 * t215 - t211, 0, 0, 0, 0; 0, qJD(2) * t202, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:10
	% EndTime: 2019-10-10 09:18:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (45->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t231 = sin(pkin(10));
	t233 = sin(qJ(1));
	t247 = t231 * t233;
	t234 = cos(qJ(1));
	t246 = t231 * t234;
	t232 = cos(pkin(10));
	t245 = t232 * t233;
	t244 = t232 * t234;
	t243 = qJD(1) * t233;
	t242 = qJD(1) * t234;
	t230 = qJ(2) + pkin(9);
	t228 = sin(t230);
	t241 = qJD(2) * t228;
	t240 = qJD(2) * t233;
	t239 = qJD(2) * t234;
	t238 = t228 * t240;
	t237 = t228 * t239;
	t229 = cos(t230);
	t236 = -t228 * t242 - t229 * t240;
	t235 = t228 * t243 - t229 * t239;
	t1 = [t232 * t238 + (-t229 * t244 - t247) * qJD(1), t235 * t232, 0, 0, 0, 0; -t232 * t237 + (-t229 * t245 + t246) * qJD(1), t236 * t232, 0, 0, 0, 0; 0, -t232 * t241, 0, 0, 0, 0; t236, -t229 * t243 - t237, 0, 0, 0, 0; -t235, t229 * t242 - t238, 0, 0, 0, 0; 0, qJD(2) * t229, 0, 0, 0, 0; t231 * t238 + (-t229 * t246 + t245) * qJD(1), t235 * t231, 0, 0, 0, 0; -t231 * t237 + (-t229 * t247 - t244) * qJD(1), t236 * t231, 0, 0, 0, 0; 0, -t231 * t241, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:11
	% EndTime: 2019-10-10 09:18:11
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (206->42), mult. (385->83), div. (0->0), fcn. (401->8), ass. (0->42)
	t327 = qJ(2) + pkin(9);
	t326 = cos(t327);
	t325 = sin(t327);
	t331 = sin(qJ(1));
	t351 = qJD(2) * t331;
	t348 = t325 * t351;
	t333 = cos(qJ(1));
	t352 = qJD(1) * t333;
	t360 = -t326 * t352 + t348;
	t328 = sin(pkin(10));
	t329 = cos(pkin(10));
	t330 = sin(qJ(6));
	t332 = cos(qJ(6));
	t343 = t328 * t330 + t329 * t332;
	t359 = qJD(6) * t343;
	t353 = qJD(1) * t331;
	t315 = -t360 * t328 - t329 * t353;
	t354 = t333 * t329;
	t357 = t331 * t328;
	t320 = t326 * t354 + t357;
	t316 = qJD(1) * t320 - t329 * t348;
	t317 = t326 * t357 + t354;
	t355 = t333 * t328;
	t356 = t331 * t329;
	t318 = t326 * t356 - t355;
	t358 = -t315 * t332 + t316 * t330 + (t317 * t330 + t318 * t332) * qJD(6);
	t350 = qJD(2) * t333;
	t347 = t325 * t350;
	t344 = t328 * t332 - t329 * t330;
	t342 = t344 * t326;
	t341 = qJD(2) * t344;
	t340 = qJD(2) * t343;
	t339 = qJD(6) * t344;
	t336 = t326 * t341;
	t335 = t326 * t340;
	t334 = -t315 * t330 - t316 * t332 + (-t317 * t332 + t318 * t330) * qJD(6);
	t319 = t326 * t355 - t356;
	t314 = -qJD(1) * t318 - t329 * t347;
	t313 = -qJD(1) * t317 - t328 * t347;
	t312 = t313 * t330 + t314 * t332 + (t319 * t332 - t320 * t330) * qJD(6);
	t311 = t313 * t332 - t314 * t330 + (-t319 * t330 - t320 * t332) * qJD(6);
	t1 = [t334, -t333 * t335 + (-t333 * t339 + t343 * t353) * t325, 0, 0, 0, t311; t312, -t331 * t335 + (-t331 * t339 - t343 * t352) * t325, 0, 0, 0, -t358; 0, qJD(6) * t342 - t325 * t340, 0, 0, 0, qJD(2) * t342 - t325 * t359; t358, -t333 * t336 + (t333 * t359 + t344 * t353) * t325, 0, 0, 0, -t312; t311, -t331 * t336 + (t331 * t359 - t344 * t352) * t325, 0, 0, 0, t334; 0, -t325 * t341 - t326 * t359, 0, 0, 0, -t325 * t339 - t335; t325 * t352 + t326 * t351, t326 * t353 + t347, 0, 0, 0, 0; t325 * t353 - t326 * t350, t360, 0, 0, 0, 0; 0, -qJD(2) * t326, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end