% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR1
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
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:59
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiRD_rot_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t31 = sin(qJ(1));
	t38 = qJD(1) * t31;
	t33 = cos(qJ(1));
	t37 = qJD(1) * t33;
	t30 = sin(qJ(2));
	t36 = qJD(2) * t30;
	t32 = cos(qJ(2));
	t35 = qJD(2) * t32;
	t34 = qJD(2) * t33;
	t29 = t31 * t36 - t32 * t37;
	t28 = t30 * t37 + t31 * t35;
	t27 = t30 * t34 + t32 * t38;
	t26 = t30 * t38 - t32 * t34;
	t1 = [t29, t26, 0, 0, 0; -t27, -t28, 0, 0, 0; 0, -t36, 0, 0, 0; t28, t27, 0, 0, 0; t26, t29, 0, 0, 0; 0, -t35, 0, 0, 0; -t38, 0, 0, 0, 0; t37, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t37 = sin(qJ(1));
	t44 = qJD(1) * t37;
	t39 = cos(qJ(1));
	t43 = qJD(1) * t39;
	t36 = sin(qJ(2));
	t42 = qJD(2) * t36;
	t38 = cos(qJ(2));
	t41 = qJD(2) * t38;
	t40 = qJD(2) * t39;
	t35 = t37 * t42 - t38 * t43;
	t34 = t36 * t43 + t37 * t41;
	t33 = t36 * t40 + t38 * t44;
	t32 = t36 * t44 - t38 * t40;
	t1 = [t35, t32, 0, 0, 0; -t33, -t34, 0, 0, 0; 0, -t42, 0, 0, 0; t34, t33, 0, 0, 0; t32, t35, 0, 0, 0; 0, -t41, 0, 0, 0; -t44, 0, 0, 0, 0; t43, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:18
	% EndTime: 2019-10-09 20:59:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t72 = qJ(2) + qJ(4);
	t69 = sin(t72);
	t71 = qJD(2) + qJD(4);
	t80 = t71 * t69;
	t70 = cos(t72);
	t79 = t71 * t70;
	t73 = sin(qJ(1));
	t78 = t71 * t73;
	t74 = cos(qJ(1));
	t77 = t71 * t74;
	t76 = qJD(1) * t73;
	t75 = qJD(1) * t74;
	t68 = t69 * t78 - t70 * t75;
	t67 = t69 * t75 + t70 * t78;
	t66 = t69 * t77 + t70 * t76;
	t65 = t69 * t76 - t70 * t77;
	t1 = [t68, t65, 0, t65, 0; -t66, -t67, 0, -t67, 0; 0, -t80, 0, -t80, 0; t67, t66, 0, t66, 0; t65, t68, 0, t68, 0; 0, -t79, 0, -t79, 0; -t76, 0, 0, 0, 0; t75, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:59:20
	% EndTime: 2019-10-09 20:59:20
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t320 = qJD(2) + qJD(4);
	t322 = sin(qJ(5));
	t344 = t320 * t322;
	t323 = sin(qJ(1));
	t343 = t320 * t323;
	t324 = cos(qJ(5));
	t342 = t320 * t324;
	t325 = cos(qJ(1));
	t341 = t320 * t325;
	t340 = t324 * t325;
	t339 = qJD(1) * t323;
	t338 = qJD(1) * t325;
	t337 = qJD(5) * t322;
	t336 = qJD(5) * t324;
	t335 = qJD(5) * t325;
	t321 = qJ(2) + qJ(4);
	t318 = sin(t321);
	t334 = t318 * t342;
	t333 = t318 * t343;
	t319 = cos(t321);
	t332 = t319 * t343;
	t331 = t318 * t341;
	t330 = t319 * t341;
	t329 = qJD(5) * t319 - qJD(1);
	t328 = qJD(1) * t319 - qJD(5);
	t327 = t329 * t322;
	t326 = t328 * t323 + t331;
	t317 = t320 * t319;
	t316 = t319 * t338 - t333;
	t315 = -t319 * t339 - t331;
	t314 = -t319 * t337 - t334;
	t313 = t318 * t344 - t319 * t336;
	t312 = -t324 * t332 + (t323 * t337 - t324 * t338) * t318;
	t311 = t322 * t332 + (t322 * t338 + t323 * t336) * t318;
	t310 = -t324 * t330 + (t322 * t335 + t324 * t339) * t318;
	t309 = t322 * t330 + (-t322 * t339 + t324 * t335) * t318;
	t308 = -t328 * t340 + (t327 + t334) * t323;
	t307 = t329 * t324 * t323 + (t328 * t325 - t333) * t322;
	t306 = t326 * t324 + t325 * t327;
	t305 = t326 * t322 - t329 * t340;
	t1 = [t308, t310, 0, t310, t305; -t306, t312, 0, t312, -t307; 0, t314, 0, t314, -t318 * t336 - t319 * t344; t307, t309, 0, t309, t306; t305, t311, 0, t311, t308; 0, t313, 0, t313, t318 * t337 - t319 * t342; -t318 * t338 - t332, t315, 0, t315, 0; -t318 * t339 + t330, t316, 0, t316, 0; 0, t317, 0, t317, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end