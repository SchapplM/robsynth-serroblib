% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:00:09
	% EndTime: 2019-12-29 18:00:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:00:09
	% EndTime: 2019-12-29 18:00:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:00:04
	% EndTime: 2019-12-29 18:00:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t11, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:00:04
	% EndTime: 2019-12-29 18:00:04
	% DurationCPUTime: 0.05s
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
	t1 = [t30, 0, t33, 0, 0; t32, 0, t31, 0, 0; 0, 0, -t39, 0, 0; -t31, 0, -t32, 0, 0; t33, 0, t30, 0, 0; 0, 0, t40, 0, 0; -t41, 0, 0, 0, 0; -t42, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:00:04
	% EndTime: 2019-12-29 18:00:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (60->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t71 = qJ(3) + qJ(4);
	t69 = cos(t71);
	t70 = qJD(3) + qJD(4);
	t78 = t70 * t69;
	t72 = sin(qJ(1));
	t77 = t70 * t72;
	t73 = cos(qJ(1));
	t76 = t70 * t73;
	t75 = qJD(1) * t72;
	t74 = qJD(1) * t73;
	t68 = sin(t71);
	t67 = t70 * t68;
	t66 = -t68 * t77 + t69 * t74;
	t65 = t68 * t74 + t69 * t77;
	t64 = t68 * t76 + t69 * t75;
	t63 = -t68 * t75 + t69 * t76;
	t1 = [t63, 0, t66, t66, 0; t65, 0, t64, t64, 0; 0, 0, -t78, -t78, 0; -t64, 0, -t65, -t65, 0; t66, 0, t63, t63, 0; 0, 0, t67, t67, 0; -t74, 0, 0, 0, 0; -t75, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:00:11
	% EndTime: 2019-12-29 18:00:12
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (166->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t328 = cos(qJ(1));
	t324 = qJ(3) + qJ(4);
	t321 = sin(t324);
	t330 = qJD(1) * t321 + qJD(5);
	t348 = t330 * t328;
	t326 = sin(qJ(1));
	t322 = cos(t324);
	t323 = qJD(3) + qJD(4);
	t342 = t323 * t328;
	t332 = t322 * t342;
	t347 = t330 * t326 - t332;
	t346 = t323 * t321;
	t325 = sin(qJ(5));
	t345 = t323 * t325;
	t344 = t323 * t326;
	t327 = cos(qJ(5));
	t343 = t323 * t327;
	t341 = qJD(1) * t326;
	t340 = qJD(1) * t328;
	t339 = qJD(5) * t325;
	t338 = qJD(5) * t327;
	t337 = qJD(5) * t328;
	t336 = t322 * t343;
	t335 = t321 * t344;
	t334 = t322 * t344;
	t333 = t321 * t342;
	t331 = -qJD(5) * t321 - qJD(1);
	t329 = t331 * t328;
	t320 = t321 * t340 + t334;
	t319 = t321 * t341 - t332;
	t318 = t321 * t339 - t336;
	t317 = t321 * t338 + t322 * t345;
	t316 = -t327 * t335 + (-t326 * t339 + t327 * t340) * t322;
	t315 = t325 * t335 + (-t325 * t340 - t326 * t338) * t322;
	t314 = t327 * t333 + (t325 * t337 + t327 * t341) * t322;
	t313 = -t325 * t333 + (-t325 * t341 + t327 * t337) * t322;
	t312 = t327 * t348 + (t331 * t325 + t336) * t326;
	t311 = t331 * t327 * t326 + (-t334 - t348) * t325;
	t310 = t325 * t329 - t347 * t327;
	t309 = t347 * t325 + t327 * t329;
	t1 = [t310, 0, t316, t316, t311; t312, 0, t314, t314, -t309; 0, 0, t318, t318, t321 * t345 - t322 * t338; t309, 0, t315, t315, -t312; t311, 0, t313, t313, t310; 0, 0, t317, t317, t321 * t343 + t322 * t339; t322 * t341 + t333, 0, t320, t320, 0; -t322 * t340 + t335, 0, t319, t319, 0; 0, 0, -t346, -t346, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end