% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiRD_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
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
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t12 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [-t11, 0, 0, 0, 0; -t12, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t12, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t33 = sin(qJ(1));
	t40 = qJD(1) * t33;
	t35 = cos(qJ(1));
	t39 = qJD(1) * t35;
	t32 = sin(qJ(3));
	t38 = qJD(3) * t32;
	t34 = cos(qJ(3));
	t37 = qJD(3) * t34;
	t36 = qJD(3) * t35;
	t31 = t33 * t38 - t34 * t39;
	t30 = t32 * t39 + t33 * t37;
	t29 = t32 * t36 + t34 * t40;
	t28 = t32 * t40 - t34 * t36;
	t1 = [t31, 0, t28, 0, 0; -t29, 0, -t30, 0, 0; 0, 0, -t38, 0, 0; t30, 0, t29, 0, 0; t28, 0, t31, 0, 0; 0, 0, -t37, 0, 0; -t40, 0, 0, 0, 0; t39, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:19
	% EndTime: 2019-12-05 18:10:19
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t226 = cos(qJ(4));
	t228 = cos(qJ(1));
	t250 = t226 * t228;
	t225 = sin(qJ(1));
	t249 = qJD(1) * t225;
	t227 = cos(qJ(3));
	t248 = qJD(1) * t227;
	t247 = qJD(1) * t228;
	t224 = sin(qJ(3));
	t246 = qJD(3) * t224;
	t245 = qJD(3) * t227;
	t244 = qJD(3) * t228;
	t223 = sin(qJ(4));
	t243 = qJD(4) * t223;
	t242 = qJD(4) * t224;
	t241 = qJD(4) * t227;
	t240 = t226 * t246;
	t239 = t226 * t242;
	t238 = t225 * t246;
	t237 = t225 * t245;
	t236 = t224 * t244;
	t235 = t227 * t244;
	t234 = -qJD(1) + t241;
	t233 = -qJD(4) + t248;
	t232 = t234 * t223;
	t231 = t224 * t247 + t237;
	t230 = -t224 * t249 + t235;
	t229 = t225 * t233 + t236;
	t222 = -t233 * t250 + (t232 + t240) * t225;
	t221 = t234 * t226 * t225 + (t228 * t233 - t238) * t223;
	t220 = t226 * t229 + t228 * t232;
	t219 = t223 * t229 - t234 * t250;
	t1 = [t222, 0, -t226 * t235 + (t226 * t249 + t228 * t243) * t224, t219, 0; -t220, 0, -t226 * t237 + (t225 * t243 - t226 * t247) * t224, -t221, 0; 0, 0, -t223 * t241 - t240, -t223 * t245 - t239, 0; t221, 0, t223 * t230 + t228 * t239, t220, 0; t219, 0, t223 * t231 + t225 * t239, t222, 0; 0, 0, t223 * t246 - t226 * t241, t223 * t242 - t226 * t245, 0; -t231, 0, -t225 * t248 - t236, 0, 0; t230, 0, t227 * t247 - t238, 0, 0; 0, 0, t245, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:20
	% EndTime: 2019-12-05 18:10:20
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (162->55), mult. (548->107), div. (0->0), fcn. (560->8), ass. (0->60)
	t373 = sin(qJ(5));
	t375 = sin(qJ(3));
	t377 = cos(qJ(5));
	t379 = cos(qJ(3));
	t378 = cos(qJ(4));
	t395 = qJD(3) * t378 - qJD(5);
	t390 = t395 * t379;
	t396 = qJD(5) * t378 - qJD(3);
	t374 = sin(qJ(4));
	t410 = qJD(4) * t374;
	t400 = t375 * t410;
	t383 = t396 * t377 * t375 + (t390 - t400) * t373;
	t408 = qJD(4) * t379;
	t402 = t374 * t408;
	t424 = t395 * t375 + t402;
	t376 = sin(qJ(1));
	t397 = qJD(1) * t379 - qJD(4);
	t380 = cos(qJ(1));
	t411 = qJD(3) * t380;
	t423 = t375 * t411 + t397 * t376;
	t412 = qJD(3) * t379;
	t414 = qJD(1) * t380;
	t387 = -t375 * t414 - t376 * t412;
	t417 = t380 * t374;
	t418 = t376 * t379;
	t385 = -qJD(5) * (t378 * t418 - t417) - t387;
	t416 = t380 * t378;
	t368 = t376 * t374 + t379 * t416;
	t409 = qJD(4) * t378;
	t399 = t380 * t409;
	t413 = qJD(3) * t375;
	t404 = t376 * t413;
	t364 = t368 * qJD(1) - t376 * t402 - t378 * t404 - t399;
	t406 = qJD(5) * t375;
	t393 = t376 * t406 + t364;
	t422 = t393 * t373 - t385 * t377;
	t420 = t375 * t378;
	t419 = t376 * t378;
	t415 = qJD(1) * t376;
	t407 = qJD(5) * t373;
	t405 = qJD(5) * t377;
	t401 = t375 * t409;
	t398 = qJD(1) - t408;
	t392 = t398 * t380;
	t362 = t374 * t392 - t423 * t378;
	t394 = t380 * t406 + t362;
	t391 = t396 * t379;
	t389 = -t373 * t379 + t377 * t420;
	t388 = t373 * t420 + t377 * t379;
	t386 = t375 * t415 - t379 * t411;
	t384 = -qJD(5) * t368 - t386;
	t382 = -t377 * t390 + (t396 * t373 + t377 * t410) * t375;
	t381 = -t385 * t373 - t393 * t377;
	t367 = -t379 * t417 + t419;
	t365 = -t374 * t418 - t416;
	t363 = t398 * t419 + (-t397 * t380 + t404) * t374;
	t361 = t423 * t374 + t378 * t392;
	t360 = t384 * t373 + t394 * t377;
	t359 = -t394 * t373 + t384 * t377;
	t1 = [t381, 0, t382 * t380 + t389 * t415, t361 * t377 - t367 * t407, t359; t360, 0, t382 * t376 - t389 * t414, t363 * t377 - t365 * t407, -t422; 0, 0, -t373 * t391 - t424 * t377, -t377 * t401 + (t373 * t406 - t377 * t412) * t374, -t383; t422, 0, t383 * t380 - t388 * t415, -t361 * t373 - t367 * t405, -t360; t359, 0, t383 * t376 + t388 * t414, -t363 * t373 - t365 * t405, t381; 0, 0, t424 * t373 - t377 * t391, t375 * t374 * t405 + (t374 * t412 + t401) * t373, t382; t363, 0, t386 * t374 - t375 * t399, t362, 0; -t361, 0, t387 * t374 - t376 * t401, t364, 0; 0, 0, -t374 * t413 + t378 * t408, t378 * t412 - t400, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end