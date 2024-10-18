% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR13
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
%   Siehe auch: S5RRRRR13_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.00s
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
	% StartTime: 2024-09-27 17:33:11
	% EndTime: 2024-09-27 17:33:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (57->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t59 = qJD(1) + qJD(2) + qJD(3);
	t60 = qJ(1) + qJ(2) + qJ(3);
	t61 = t59 * cos(t60);
	t56 = t59 * sin(t60);
	t1 = [-t61, -t61, -t61, 0, 0; -t56, -t56, -t56, 0, 0; 0, 0, 0, 0, 0; t56, t56, t56, 0, 0; -t61, -t61, -t61, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:12
	% EndTime: 2024-09-27 17:33:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (269->19), mult. (176->24), div. (0->0), fcn. (176->6), ass. (0->23)
	t205 = qJD(1) + qJD(2) + qJD(3);
	t207 = sin(pkin(5));
	t219 = t205 * t207;
	t208 = cos(pkin(5));
	t209 = sin(qJ(4));
	t218 = t208 * t209;
	t210 = cos(qJ(4));
	t217 = t208 * t210;
	t216 = qJD(4) * t207;
	t206 = qJ(1) + qJ(2) + qJ(3);
	t203 = sin(t206);
	t215 = t203 * t219;
	t204 = cos(t206);
	t214 = t203 * t209 - t204 * t217;
	t213 = t203 * t210 + t204 * t218;
	t212 = t203 * t217 + t204 * t209;
	t211 = t203 * t218 - t204 * t210;
	t202 = t204 * t219;
	t201 = t214 * qJD(4) + t211 * t205;
	t200 = t213 * qJD(4) + t212 * t205;
	t199 = t212 * qJD(4) + t213 * t205;
	t198 = t211 * qJD(4) + t214 * t205;
	t1 = [t201, t201, t201, t198, 0; -t199, -t199, -t199, -t200, 0; 0, 0, 0, -t209 * t216, 0; t200, t200, t200, t199, 0; t198, t198, t198, t201, 0; 0, 0, 0, -t210 * t216, 0; -t215, -t215, -t215, 0, 0; t202, t202, t202, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 17:33:12
	% EndTime: 2024-09-27 17:33:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (655->27), mult. (288->27), div. (0->0), fcn. (220->9), ass. (0->32)
	t332 = qJD(4) + qJD(5);
	t348 = -t332 / 0.2e1;
	t333 = qJ(4) + qJ(5);
	t328 = pkin(5) - t333;
	t322 = sin(t328);
	t352 = t322 * t348;
	t326 = qJD(1) + qJD(2) + qJD(3);
	t329 = sin(t333);
	t327 = pkin(5) + t333;
	t346 = sin(t327);
	t351 = t329 * t332 + (t346 / 0.2e1 - t322 / 0.2e1) * t326;
	t317 = t346 * t348;
	t350 = t326 * t329 - t317 + t352;
	t331 = qJ(1) + qJ(2) + qJ(3);
	t324 = sin(t331);
	t325 = cos(t331);
	t347 = cos(t328);
	t349 = cos(t327) / 0.2e1;
	t315 = t347 / 0.2e1 + t349;
	t330 = cos(t333);
	t339 = t315 * t332 + t326 * t330;
	t304 = t339 * t324 + t351 * t325;
	t343 = t326 * sin(pkin(5));
	t340 = t324 * t343;
	t336 = t326 * t315 + t332 * t330;
	t335 = t351 * t324 - t339 * t325;
	t316 = t325 * t343;
	t311 = t332 * t349 + t347 * t348;
	t310 = t317 + t352;
	t306 = t336 * t324 + t350 * t325;
	t303 = t350 * t324 - t336 * t325;
	t1 = [t335, t335, t335, t303, t303; -t304, -t304, -t304, -t306, -t306; 0, 0, 0, t311, t311; t306, t306, t306, t304, t304; t303, t303, t303, t335, t335; 0, 0, 0, t310, t310; -t340, -t340, -t340, 0, 0; t316, t316, t316, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end