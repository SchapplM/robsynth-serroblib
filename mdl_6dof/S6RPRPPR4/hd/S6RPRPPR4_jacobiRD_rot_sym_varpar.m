% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
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
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
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
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
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
	% StartTime: 2019-10-10 00:20:30
	% EndTime: 2019-10-10 00:20:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->18), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t194 = sin(pkin(10));
	t196 = sin(qJ(1));
	t210 = t194 * t196;
	t197 = cos(qJ(1));
	t209 = t194 * t197;
	t195 = cos(pkin(10));
	t208 = t195 * t196;
	t207 = t195 * t197;
	t206 = qJD(1) * t196;
	t205 = qJD(1) * t197;
	t193 = pkin(9) + qJ(3);
	t191 = sin(t193);
	t204 = qJD(3) * t191;
	t203 = qJD(3) * t196;
	t202 = qJD(3) * t197;
	t201 = t191 * t203;
	t200 = t191 * t202;
	t192 = cos(t193);
	t199 = t191 * t205 + t192 * t203;
	t198 = t191 * t206 - t192 * t202;
	t1 = [t195 * t201 + (-t192 * t207 - t210) * qJD(1), 0, t198 * t195, 0, 0, 0; -t195 * t200 + (-t192 * t208 + t209) * qJD(1), 0, -t199 * t195, 0, 0, 0; 0, 0, -t195 * t204, 0, 0, 0; -t194 * t201 + (t192 * t209 - t208) * qJD(1), 0, -t198 * t194, 0, 0, 0; t194 * t200 + (t192 * t210 + t207) * qJD(1), 0, t199 * t194, 0, 0, 0; 0, 0, t194 * t204, 0, 0, 0; -t199, 0, -t192 * t206 - t200, 0, 0, 0; -t198, 0, t192 * t205 - t201, 0, 0, 0; 0, 0, qJD(3) * t192, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:30
	% EndTime: 2019-10-10 00:20:30
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (45->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t227 = sin(pkin(10));
	t229 = sin(qJ(1));
	t243 = t227 * t229;
	t230 = cos(qJ(1));
	t242 = t227 * t230;
	t228 = cos(pkin(10));
	t241 = t228 * t229;
	t240 = t228 * t230;
	t239 = qJD(1) * t229;
	t238 = qJD(1) * t230;
	t226 = pkin(9) + qJ(3);
	t224 = sin(t226);
	t237 = qJD(3) * t224;
	t236 = qJD(3) * t229;
	t235 = qJD(3) * t230;
	t234 = t224 * t236;
	t233 = t224 * t235;
	t225 = cos(t226);
	t232 = -t224 * t238 - t225 * t236;
	t231 = t224 * t239 - t225 * t235;
	t1 = [t228 * t234 + (-t225 * t240 - t243) * qJD(1), 0, t231 * t228, 0, 0, 0; -t228 * t233 + (-t225 * t241 + t242) * qJD(1), 0, t232 * t228, 0, 0, 0; 0, 0, -t228 * t237, 0, 0, 0; t232, 0, -t225 * t239 - t233, 0, 0, 0; -t231, 0, t225 * t238 - t234, 0, 0, 0; 0, 0, qJD(3) * t225, 0, 0, 0; t227 * t234 + (-t225 * t242 + t241) * qJD(1), 0, t231 * t227, 0, 0, 0; -t227 * t233 + (-t225 * t243 - t240) * qJD(1), 0, t232 * t227, 0, 0, 0; 0, 0, -t227 * t237, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:30
	% EndTime: 2019-10-10 00:20:30
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (206->42), mult. (385->83), div. (0->0), fcn. (401->8), ass. (0->42)
	t319 = pkin(9) + qJ(3);
	t318 = cos(t319);
	t317 = sin(t319);
	t323 = sin(qJ(1));
	t343 = qJD(3) * t323;
	t340 = t317 * t343;
	t325 = cos(qJ(1));
	t344 = qJD(1) * t325;
	t352 = -t318 * t344 + t340;
	t320 = sin(pkin(10));
	t321 = cos(pkin(10));
	t322 = sin(qJ(6));
	t324 = cos(qJ(6));
	t335 = t320 * t322 + t321 * t324;
	t351 = qJD(6) * t335;
	t345 = qJD(1) * t323;
	t307 = -t352 * t320 - t321 * t345;
	t346 = t325 * t321;
	t349 = t323 * t320;
	t312 = t318 * t346 + t349;
	t308 = t312 * qJD(1) - t321 * t340;
	t309 = t318 * t349 + t346;
	t347 = t325 * t320;
	t348 = t323 * t321;
	t310 = t318 * t348 - t347;
	t350 = -t307 * t324 + t308 * t322 + (t309 * t322 + t310 * t324) * qJD(6);
	t342 = qJD(3) * t325;
	t339 = t317 * t342;
	t336 = t320 * t324 - t321 * t322;
	t334 = t336 * t318;
	t333 = qJD(3) * t336;
	t332 = qJD(3) * t335;
	t331 = qJD(6) * t336;
	t328 = t318 * t333;
	t327 = t318 * t332;
	t326 = -t307 * t322 - t308 * t324 + (-t309 * t324 + t310 * t322) * qJD(6);
	t311 = t318 * t347 - t348;
	t306 = -t310 * qJD(1) - t321 * t339;
	t305 = -t309 * qJD(1) - t320 * t339;
	t304 = t305 * t322 + t306 * t324 + (t311 * t324 - t312 * t322) * qJD(6);
	t303 = t305 * t324 - t306 * t322 + (-t311 * t322 - t312 * t324) * qJD(6);
	t1 = [t326, 0, -t325 * t327 + (-t325 * t331 + t335 * t345) * t317, 0, 0, t303; t304, 0, -t323 * t327 + (-t323 * t331 - t335 * t344) * t317, 0, 0, -t350; 0, 0, qJD(6) * t334 - t317 * t332, 0, 0, qJD(3) * t334 - t317 * t351; t350, 0, -t325 * t328 + (t325 * t351 + t336 * t345) * t317, 0, 0, -t304; t303, 0, -t323 * t328 + (t323 * t351 - t336 * t344) * t317, 0, 0, t326; 0, 0, -t317 * t333 - t318 * t351, 0, 0, -t317 * t331 - t327; t317 * t344 + t318 * t343, 0, t318 * t345 + t339, 0, 0, 0; t317 * t345 - t318 * t342, 0, t352, 0, 0, 0; 0, 0, -qJD(3) * t318, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end