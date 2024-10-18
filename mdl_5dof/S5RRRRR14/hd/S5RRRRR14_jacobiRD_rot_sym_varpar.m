% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRRRR14_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR14_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR14_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:41
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:41
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0, 0; -t44, -t44, 0, 0, 0; 0, 0, 0, 0, 0; t44, t44, 0, 0, 0; -t49, -t49, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (120->17), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t193 = qJD(1) + qJD(2);
	t195 = sin(pkin(5));
	t207 = t193 * t195;
	t196 = cos(pkin(5));
	t197 = sin(qJ(3));
	t206 = t196 * t197;
	t198 = cos(qJ(3));
	t205 = t196 * t198;
	t204 = qJD(3) * t195;
	t194 = qJ(1) + qJ(2);
	t191 = sin(t194);
	t203 = t191 * t207;
	t192 = cos(t194);
	t202 = t191 * t197 - t192 * t205;
	t201 = t191 * t198 + t192 * t206;
	t200 = t191 * t205 + t192 * t197;
	t199 = t191 * t206 - t192 * t198;
	t190 = t192 * t207;
	t189 = t202 * qJD(3) + t199 * t193;
	t188 = t201 * qJD(3) + t200 * t193;
	t187 = t200 * qJD(3) + t201 * t193;
	t186 = t199 * qJD(3) + t202 * t193;
	t1 = [t189, t189, t186, 0, 0; -t187, -t187, -t188, 0, 0; 0, 0, -t197 * t204, 0, 0; t188, t188, t187, 0, 0; t186, t186, t189, 0, 0; 0, 0, -t198 * t204, 0, 0; -t203, -t203, 0, 0, 0; t190, t190, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (422->25), mult. (232->27), div. (0->0), fcn. (176->9), ass. (0->32)
	t316 = qJD(3) + qJD(4);
	t334 = -t316 / 0.2e1;
	t318 = qJ(3) + qJ(4);
	t311 = pkin(5) - t318;
	t308 = sin(t311);
	t338 = t308 * t334;
	t312 = sin(t318);
	t317 = qJD(1) + qJD(2);
	t310 = pkin(5) + t318;
	t332 = sin(t310);
	t337 = t312 * t316 + (t332 / 0.2e1 - t308 / 0.2e1) * t317;
	t302 = t332 * t334;
	t336 = t317 * t312 - t302 + t338;
	t319 = qJ(1) + qJ(2);
	t313 = sin(t319);
	t315 = cos(t319);
	t333 = cos(t311);
	t335 = cos(t310) / 0.2e1;
	t300 = t333 / 0.2e1 + t335;
	t314 = cos(t318);
	t325 = t300 * t316 + t314 * t317;
	t290 = t325 * t313 + t337 * t315;
	t327 = t317 * sin(pkin(5));
	t326 = t313 * t327;
	t322 = t317 * t300 + t316 * t314;
	t321 = t337 * t313 - t325 * t315;
	t305 = t315 * t327;
	t297 = t316 * t335 + t333 * t334;
	t296 = t302 + t338;
	t292 = t322 * t313 + t336 * t315;
	t289 = t336 * t313 - t322 * t315;
	t1 = [t321, t321, t289, t289, 0; -t290, -t290, -t292, -t292, 0; 0, 0, t297, t297, 0; t292, t292, t290, t290, 0; t289, t289, t321, t321, 0; 0, 0, t296, t296, 0; -t326, -t326, 0, 0, 0; t305, t305, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:44:42
	% EndTime: 2024-09-27 18:44:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (736->25), mult. (292->27), div. (0->0), fcn. (220->9), ass. (0->30)
	t351 = qJ(3) + qJ(4) + qJ(5);
	t346 = sin(t351);
	t348 = qJD(3) + qJD(4) + qJD(5);
	t352 = qJD(1) + qJD(2);
	t344 = pkin(5) + t351;
	t361 = -sin(t344) / 0.2e1;
	t345 = pkin(5) - t351;
	t368 = sin(t345);
	t357 = t368 / 0.2e1 + t361;
	t373 = t346 * t348 - t357 * t352;
	t372 = t352 * t346 - t357 * t348;
	t369 = cos(t345);
	t370 = cos(t344) / 0.2e1;
	t336 = t369 / 0.2e1 + t370;
	t347 = cos(t351);
	t371 = t352 * t336 + t348 * t347;
	t353 = qJ(1) + qJ(2);
	t349 = sin(t353);
	t350 = cos(t353);
	t359 = t336 * t348 + t347 * t352;
	t326 = t359 * t349 + t373 * t350;
	t362 = t352 * sin(pkin(5));
	t360 = t349 * t362;
	t355 = t373 * t349 - t359 * t350;
	t340 = t350 * t362;
	t333 = (t370 - t369 / 0.2e1) * t348;
	t332 = (t361 - t368 / 0.2e1) * t348;
	t328 = t371 * t349 + t372 * t350;
	t325 = t372 * t349 - t371 * t350;
	t1 = [t355, t355, t325, t325, t325; -t326, -t326, -t328, -t328, -t328; 0, 0, t333, t333, t333; t328, t328, t326, t326, t326; t325, t325, t355, t355, t355; 0, 0, t332, t332, t332; -t360, -t360, 0, 0, 0; t340, t340, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end