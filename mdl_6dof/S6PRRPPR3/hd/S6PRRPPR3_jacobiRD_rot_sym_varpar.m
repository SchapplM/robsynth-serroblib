% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(10));
	t58 = sin(pkin(10));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:55
	% EndTime: 2019-10-09 22:10:55
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(6));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(6));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(10));
	t212 = cos(pkin(10));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0, 0; 0, t204, 0, 0, 0, 0; 0, t202, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:56
	% EndTime: 2019-10-09 22:10:56
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t269 = sin(pkin(6));
	t272 = sin(qJ(3));
	t285 = t269 * t272;
	t274 = cos(qJ(3));
	t284 = t269 * t274;
	t271 = cos(pkin(6));
	t273 = sin(qJ(2));
	t283 = t271 * t273;
	t275 = cos(qJ(2));
	t282 = t271 * t275;
	t281 = qJD(2) * t273;
	t280 = qJD(3) * t272;
	t279 = qJD(3) * t274;
	t278 = qJD(3) * t275;
	t277 = t269 * qJD(2) * t275;
	t268 = sin(pkin(10));
	t270 = cos(pkin(10));
	t264 = -t268 * t273 + t270 * t282;
	t265 = t268 * t275 + t270 * t283;
	t266 = -t268 * t282 - t270 * t273;
	t276 = t268 * t283 - t270 * t275;
	t263 = t276 * qJD(2);
	t262 = t266 * qJD(2);
	t261 = t265 * qJD(2);
	t260 = t264 * qJD(2);
	t1 = [0, t263 * t274 - t266 * t280, -t262 * t272 + (-t268 * t285 + t274 * t276) * qJD(3), 0, 0, 0; 0, -t261 * t274 - t264 * t280, -t260 * t272 + (-t265 * t274 + t270 * t285) * qJD(3), 0, 0, 0; 0, (-t272 * t278 - t274 * t281) * t269, -t272 * t277 + (-t271 * t272 - t273 * t284) * qJD(3), 0, 0, 0; 0, t262, 0, 0, 0, 0; 0, t260, 0, 0, 0, 0; 0, t277, 0, 0, 0, 0; 0, t263 * t272 + t266 * t279, t262 * t274 + (t268 * t284 + t272 * t276) * qJD(3), 0, 0, 0; 0, -t261 * t272 + t264 * t279, t260 * t274 + (-t265 * t272 - t270 * t284) * qJD(3), 0, 0, 0; 0, (-t272 * t281 + t274 * t278) * t269, t274 * t277 + (t271 * t274 - t273 * t285) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:55
	% EndTime: 2019-10-09 22:10:55
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (37->24), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->27)
	t234 = sin(pkin(6));
	t237 = sin(qJ(3));
	t252 = t234 * t237;
	t239 = cos(qJ(3));
	t251 = t234 * t239;
	t236 = cos(pkin(6));
	t238 = sin(qJ(2));
	t250 = t236 * t238;
	t240 = cos(qJ(2));
	t249 = t236 * t240;
	t248 = qJD(2) * t238;
	t247 = qJD(3) * t237;
	t246 = qJD(3) * t239;
	t245 = qJD(3) * t240;
	t235 = cos(pkin(10));
	t244 = t235 * t249;
	t243 = t234 * qJD(2) * t240;
	t233 = sin(pkin(10));
	t229 = t233 * t240 + t235 * t250;
	t242 = t233 * t249 + t235 * t238;
	t241 = t233 * t250 - t235 * t240;
	t228 = -t233 * t238 + t244;
	t227 = t241 * qJD(2);
	t226 = t242 * qJD(2);
	t225 = t229 * qJD(2);
	t224 = -qJD(2) * t244 + t233 * t248;
	t1 = [0, t227 * t237 - t242 * t246, -t226 * t239 + (t233 * t251 + t237 * t241) * qJD(3), 0, 0, 0; 0, -t225 * t237 + t228 * t246, -t224 * t239 + (-t229 * t237 - t235 * t251) * qJD(3), 0, 0, 0; 0, (-t237 * t248 + t239 * t245) * t234, t239 * t243 + (t236 * t239 - t238 * t252) * qJD(3), 0, 0, 0; 0, -t227 * t239 - t242 * t247, -t226 * t237 + (t233 * t252 - t239 * t241) * qJD(3), 0, 0, 0; 0, t225 * t239 + t228 * t247, -t224 * t237 + (t229 * t239 - t235 * t252) * qJD(3), 0, 0, 0; 0, (t237 * t245 + t239 * t248) * t234, t237 * t243 + (t236 * t237 + t238 * t251) * qJD(3), 0, 0, 0; 0, t226, 0, 0, 0, 0; 0, t224, 0, 0, 0, 0; 0, -t243, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:56
	% EndTime: 2019-10-09 22:10:57
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (153->62), mult. (516->126), div. (0->0), fcn. (564->10), ass. (0->52)
	t359 = sin(qJ(3));
	t360 = sin(qJ(2));
	t362 = cos(qJ(3));
	t363 = cos(qJ(2));
	t380 = qJD(3) * t363;
	t389 = (qJD(2) * t359 + qJD(6)) * t360 - t362 * t380;
	t355 = sin(pkin(6));
	t388 = t355 * t359;
	t387 = t355 * t362;
	t386 = t355 * t363;
	t357 = cos(pkin(6));
	t385 = t357 * t360;
	t384 = t357 * t363;
	t383 = qJD(2) * t360;
	t382 = qJD(3) * t359;
	t381 = qJD(3) * t362;
	t358 = sin(qJ(6));
	t379 = qJD(6) * t358;
	t378 = qJD(6) * t359;
	t361 = cos(qJ(6));
	t377 = qJD(6) * t361;
	t356 = cos(pkin(10));
	t376 = t356 * t384;
	t375 = t355 * t383;
	t374 = qJD(2) * t386;
	t354 = sin(pkin(10));
	t342 = -qJD(2) * t376 + t354 * t383;
	t346 = -t354 * t360 + t376;
	t371 = -t346 * t378 + t342;
	t367 = t354 * t384 + t356 * t360;
	t344 = t367 * qJD(2);
	t370 = t367 * t378 + t344;
	t369 = (-qJD(2) - t378) * t363;
	t347 = t354 * t363 + t356 * t385;
	t336 = t347 * t359 + t356 * t387;
	t337 = t347 * t362 - t356 * t388;
	t366 = t354 * t385 - t356 * t363;
	t368 = t354 * t387 + t359 * t366;
	t339 = t354 * t388 - t362 * t366;
	t351 = t357 * t359 + t360 * t387;
	t350 = -t357 * t362 + t360 * t388;
	t343 = t347 * qJD(2);
	t365 = -qJD(6) * t347 - t343 * t359 + t346 * t381;
	t345 = t366 * qJD(2);
	t364 = qJD(6) * t366 + t345 * t359 - t367 * t381;
	t341 = -qJD(3) * t350 + t362 * t374;
	t340 = qJD(3) * t351 + t359 * t374;
	t335 = qJD(3) * t368 - t344 * t362;
	t334 = qJD(3) * t339 - t344 * t359;
	t333 = -qJD(3) * t336 - t342 * t362;
	t332 = qJD(3) * t337 - t342 * t359;
	t1 = [0, t358 * t370 + t364 * t361, t335 * t361 - t339 * t379, 0, 0, -t334 * t358 + t345 * t361 + (t358 * t367 + t361 * t368) * qJD(6); 0, t371 * t358 + t365 * t361, t333 * t361 - t337 * t379, 0, 0, -t332 * t358 - t343 * t361 + (-t336 * t361 - t346 * t358) * qJD(6); 0, (t358 * t369 - t389 * t361) * t355, t341 * t361 - t351 * t379, 0, 0, -t361 * t375 - t340 * t358 + (-t350 * t361 - t358 * t386) * qJD(6); 0, -t364 * t358 + t361 * t370, -t335 * t358 - t339 * t377, 0, 0, -t334 * t361 - t345 * t358 + (-t358 * t368 + t361 * t367) * qJD(6); 0, -t365 * t358 + t371 * t361, -t333 * t358 - t337 * t377, 0, 0, -t332 * t361 + t343 * t358 + (t336 * t358 - t346 * t361) * qJD(6); 0, (t389 * t358 + t361 * t369) * t355, -t341 * t358 - t351 * t377, 0, 0, t358 * t375 - t340 * t361 + (t350 * t358 - t361 * t386) * qJD(6); 0, t345 * t362 + t367 * t382, -t334, 0, 0, 0; 0, -t343 * t362 - t346 * t382, -t332, 0, 0, 0; 0, (-t359 * t380 - t362 * t383) * t355, -t340, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end