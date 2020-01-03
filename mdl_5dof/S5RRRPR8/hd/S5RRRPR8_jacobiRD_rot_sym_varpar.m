% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:20
	% EndTime: 2019-12-31 21:21:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:21
	% EndTime: 2019-12-31 21:21:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:21
	% EndTime: 2019-12-31 21:21:21
	% DurationCPUTime: 0.03s
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
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:21
	% EndTime: 2019-12-31 21:21:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t71 = qJ(2) + qJ(3);
	t68 = sin(t71);
	t70 = qJD(2) + qJD(3);
	t79 = t70 * t68;
	t69 = cos(t71);
	t78 = t70 * t69;
	t72 = sin(qJ(1));
	t77 = t70 * t72;
	t73 = cos(qJ(1));
	t76 = t70 * t73;
	t75 = qJD(1) * t72;
	t74 = qJD(1) * t73;
	t67 = t68 * t77 - t69 * t74;
	t66 = t68 * t74 + t69 * t77;
	t65 = t68 * t76 + t69 * t75;
	t64 = t68 * t75 - t69 * t76;
	t1 = [t67, t64, t64, 0, 0; -t65, -t66, -t66, 0, 0; 0, -t79, -t79, 0, 0; t66, t65, t65, 0, 0; t64, t67, t67, 0, 0; 0, -t78, -t78, 0, 0; -t75, 0, 0, 0, 0; t74, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:21
	% EndTime: 2019-12-31 21:21:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t235 = qJD(2) + qJD(3);
	t237 = sin(qJ(1));
	t242 = t235 * t237;
	t238 = cos(qJ(1));
	t241 = t235 * t238;
	t240 = qJD(1) * t237;
	t239 = qJD(1) * t238;
	t236 = qJ(2) + qJ(3);
	t234 = cos(t236);
	t233 = sin(t236);
	t232 = t235 * t234;
	t231 = t235 * t233;
	t230 = -t233 * t242 + t234 * t239;
	t229 = t233 * t239 + t234 * t242;
	t228 = t233 * t241 + t234 * t240;
	t227 = -t233 * t240 + t234 * t241;
	t1 = [-t240, 0, 0, 0, 0; t239, 0, 0, 0, 0; 0, 0, 0, 0, 0; t230, t227, t227, 0, 0; t228, t229, t229, 0, 0; 0, t231, t231, 0, 0; -t229, -t228, -t228, 0, 0; t227, t230, t230, 0, 0; 0, t232, t232, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:22
	% EndTime: 2019-12-31 21:21:22
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (166->31), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t334 = cos(qJ(5));
	t331 = qJ(2) + qJ(3);
	t328 = sin(t331);
	t339 = qJD(5) * t328 + qJD(1);
	t356 = t334 * t339;
	t332 = sin(qJ(5));
	t355 = t339 * t332;
	t330 = qJD(2) + qJD(3);
	t354 = t330 * t328;
	t353 = t330 * t332;
	t333 = sin(qJ(1));
	t352 = t330 * t333;
	t351 = t330 * t334;
	t335 = cos(qJ(1));
	t350 = t330 * t335;
	t349 = qJD(1) * t333;
	t348 = qJD(1) * t335;
	t347 = qJD(5) * t332;
	t346 = qJD(5) * t334;
	t345 = qJD(5) * t335;
	t329 = cos(t331);
	t344 = t329 * t351;
	t343 = t328 * t352;
	t342 = t329 * t352;
	t341 = t328 * t350;
	t340 = t329 * t350;
	t338 = -qJD(1) * t328 - qJD(5);
	t337 = t338 * t335;
	t336 = t338 * t333 + t340;
	t327 = -t328 * t348 - t342;
	t326 = t328 * t349 - t340;
	t325 = -t328 * t347 + t344;
	t324 = t328 * t346 + t329 * t353;
	t323 = -t334 * t343 + (-t333 * t347 + t334 * t348) * t329;
	t322 = -t332 * t343 + (t332 * t348 + t333 * t346) * t329;
	t321 = -t334 * t341 + (-t332 * t345 - t334 * t349) * t329;
	t320 = -t332 * t341 + (-t332 * t349 + t334 * t345) * t329;
	t319 = t336 * t332 + t335 * t356;
	t318 = t336 * t334 - t335 * t355;
	t317 = -t333 * t356 + (t337 - t342) * t332;
	t316 = t334 * t337 + (-t344 + t355) * t333;
	t1 = [t317, t320, t320, 0, t318; t319, t322, t322, 0, -t316; 0, t324, t324, 0, t328 * t351 + t329 * t347; t316, t321, t321, 0, -t319; t318, t323, t323, 0, t317; 0, t325, t325, 0, -t328 * t353 + t329 * t346; -t329 * t348 + t343, t326, t326, 0, 0; -t329 * t349 - t341, t327, t327, 0, 0; 0, -t354, -t354, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end