% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:40
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t139 = cos(pkin(6));
	t140 = sin(qJ(2));
	t144 = t139 * t140;
	t141 = cos(qJ(2));
	t143 = t139 * t141;
	t142 = qJD(2) * sin(pkin(6));
	t138 = cos(pkin(10));
	t136 = sin(pkin(10));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, (-t136 * t144 + t138 * t141) * qJD(2), 0, 0, 0, 0; 0, (t136 * t141 + t138 * t144) * qJD(2), 0, 0, 0, 0; 0, t140 * t142, 0, 0, 0, 0; 0, (-t136 * t143 - t138 * t140) * qJD(2), 0, 0, 0, 0; 0, (-t136 * t140 + t138 * t143) * qJD(2), 0, 0, 0, 0; 0, t141 * t142, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:57
	% EndTime: 2019-10-09 21:40:57
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (37->25), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t209 = sin(pkin(6));
	t212 = sin(qJ(4));
	t225 = t209 * t212;
	t214 = cos(qJ(4));
	t224 = t209 * t214;
	t215 = cos(qJ(2));
	t223 = t209 * t215;
	t211 = cos(pkin(6));
	t213 = sin(qJ(2));
	t222 = t211 * t213;
	t221 = t211 * t215;
	t220 = qJD(2) * t215;
	t219 = qJD(4) * t212;
	t218 = qJD(4) * t214;
	t217 = t209 * qJD(2) * t213;
	t208 = sin(pkin(10));
	t210 = cos(pkin(10));
	t216 = -t208 * t213 + t210 * t221;
	t205 = t208 * t215 + t210 * t222;
	t206 = t208 * t221 + t210 * t213;
	t207 = -t208 * t222 + t210 * t215;
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t201 = t205 * qJD(2);
	t200 = t216 * qJD(2);
	t1 = [0, -t202 * t212 + t207 * t218, 0, t203 * t214 + (-t206 * t212 - t208 * t224) * qJD(4), 0, 0; 0, t200 * t212 + t205 * t218, 0, t201 * t214 + (t210 * t224 + t212 * t216) * qJD(4), 0, 0; 0, (t212 * t220 + t213 * t218) * t209, 0, t214 * t217 + (-t211 * t214 + t212 * t223) * qJD(4), 0, 0; 0, -t202 * t214 - t207 * t219, 0, -t203 * t212 + (-t206 * t214 + t208 * t225) * qJD(4), 0, 0; 0, t200 * t214 - t205 * t219, 0, -t201 * t212 + (-t210 * t225 + t214 * t216) * qJD(4), 0, 0; 0, (-t213 * t219 + t214 * t220) * t209, 0, -t212 * t217 + (t211 * t212 + t214 * t223) * qJD(4), 0, 0; 0, -t203, 0, 0, 0, 0; 0, -t201, 0, 0, 0, 0; 0, -t217, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:57
	% EndTime: 2019-10-09 21:40:58
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (37->25), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t259 = sin(pkin(6));
	t262 = sin(qJ(4));
	t275 = t259 * t262;
	t264 = cos(qJ(4));
	t274 = t259 * t264;
	t265 = cos(qJ(2));
	t273 = t259 * t265;
	t261 = cos(pkin(6));
	t263 = sin(qJ(2));
	t272 = t261 * t263;
	t271 = t261 * t265;
	t270 = qJD(2) * t265;
	t269 = qJD(4) * t262;
	t268 = qJD(4) * t264;
	t267 = t259 * qJD(2) * t263;
	t258 = sin(pkin(10));
	t260 = cos(pkin(10));
	t266 = -t258 * t263 + t260 * t271;
	t255 = t258 * t265 + t260 * t272;
	t256 = t258 * t271 + t260 * t263;
	t257 = -t258 * t272 + t260 * t265;
	t253 = t257 * qJD(2);
	t252 = t256 * qJD(2);
	t251 = t255 * qJD(2);
	t250 = t266 * qJD(2);
	t1 = [0, -t253, 0, 0, 0, 0; 0, -t251, 0, 0, 0, 0; 0, -t267, 0, 0, 0, 0; 0, t252 * t262 - t257 * t268, 0, -t253 * t264 + (t256 * t262 + t258 * t274) * qJD(4), 0, 0; 0, -t250 * t262 - t255 * t268, 0, -t251 * t264 + (-t260 * t274 - t262 * t266) * qJD(4), 0, 0; 0, (-t262 * t270 - t263 * t268) * t259, 0, -t264 * t267 + (t261 * t264 - t262 * t273) * qJD(4), 0, 0; 0, t252 * t264 + t257 * t269, 0, t253 * t262 + (t256 * t264 - t258 * t275) * qJD(4), 0, 0; 0, -t250 * t264 + t255 * t269, 0, t251 * t262 + (t260 * t275 - t264 * t266) * qJD(4), 0, 0; 0, (t263 * t269 - t264 * t270) * t259, 0, t262 * t267 + (-t261 * t262 - t264 * t273) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:58
	% EndTime: 2019-10-09 21:40:58
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (153->63), mult. (516->130), div. (0->0), fcn. (564->10), ass. (0->53)
	t357 = sin(pkin(6));
	t361 = sin(qJ(4));
	t389 = t357 * t361;
	t362 = sin(qJ(2));
	t388 = t357 * t362;
	t364 = cos(qJ(4));
	t387 = t357 * t364;
	t365 = cos(qJ(2));
	t386 = t357 * t365;
	t359 = cos(pkin(6));
	t385 = t359 * t362;
	t384 = t359 * t365;
	t363 = cos(qJ(6));
	t383 = t362 * t363;
	t382 = qJD(2) * t365;
	t381 = qJD(4) * t361;
	t380 = qJD(4) * t364;
	t360 = sin(qJ(6));
	t379 = qJD(6) * t360;
	t378 = qJD(6) * t363;
	t377 = qJD(6) * t364;
	t356 = sin(pkin(10));
	t376 = t356 * t385;
	t375 = qJD(2) * t388;
	t374 = t357 * t382;
	t373 = qJD(2) + t377;
	t358 = cos(pkin(10));
	t349 = t356 * t365 + t358 * t385;
	t345 = t349 * qJD(2);
	t372 = t349 * t377 + t345;
	t347 = -qJD(2) * t376 + t358 * t382;
	t351 = t358 * t365 - t376;
	t371 = t351 * t377 + t347;
	t370 = (-qJD(2) * t364 - qJD(6)) * t365;
	t368 = -t356 * t362 + t358 * t384;
	t369 = t358 * t389 - t364 * t368;
	t350 = t356 * t384 + t358 * t362;
	t338 = t350 * t361 + t356 * t387;
	t337 = -t350 * t364 + t356 * t389;
	t352 = t359 * t361 + t364 * t386;
	t353 = t359 * t364 - t361 * t386;
	t344 = t368 * qJD(2);
	t367 = -qJD(6) * t368 - t344 * t364 + t349 * t381;
	t346 = t350 * qJD(2);
	t366 = qJD(6) * t350 + t346 * t364 + t351 * t381;
	t342 = qJD(4) * t353 - t364 * t375;
	t341 = -qJD(4) * t352 + t361 * t375;
	t340 = -t358 * t387 - t361 * t368;
	t336 = t368 * t381 + (qJD(4) * t357 * t358 + t345) * t364;
	t335 = qJD(4) * t369 + t345 * t361;
	t334 = qJD(4) * t338 - t347 * t364;
	t333 = -qJD(4) * t337 + t347 * t361;
	t1 = [0, t366 * t360 - t363 * t371, 0, t333 * t360 + t338 * t378, 0, t334 * t363 + t346 * t360 + (-t337 * t360 - t351 * t363) * qJD(6); 0, t367 * t360 - t372 * t363, 0, t335 * t360 + t340 * t378, 0, -t336 * t363 - t344 * t360 + (-t349 * t363 + t360 * t369) * qJD(6); 0, (-t373 * t383 + (t362 * t381 + t370) * t360) * t357, 0, t341 * t360 + t353 * t378, 0, -t360 * t374 + t342 * t363 + (-t352 * t360 - t357 * t383) * qJD(6); 0, t360 * t371 + t366 * t363, 0, t333 * t363 - t338 * t379, 0, -t334 * t360 + t346 * t363 + (-t337 * t363 + t351 * t360) * qJD(6); 0, t360 * t372 + t367 * t363, 0, t335 * t363 - t340 * t379, 0, t336 * t360 - t344 * t363 + (t349 * t360 + t363 * t369) * qJD(6); 0, (t363 * t370 + (t373 * t360 + t363 * t381) * t362) * t357, 0, t341 * t363 - t353 * t379, 0, -t363 * t374 - t342 * t360 + (-t352 * t363 + t360 * t388) * qJD(6); 0, -t346 * t361 + t351 * t380, 0, -t334, 0, 0; 0, t344 * t361 + t349 * t380, 0, t336, 0, 0; 0, (t361 * t382 + t362 * t380) * t357, 0, -t342, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end