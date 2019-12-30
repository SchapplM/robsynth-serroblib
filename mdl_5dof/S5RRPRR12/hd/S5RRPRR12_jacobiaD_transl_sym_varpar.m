% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR12
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:19:21
	% EndTime: 2019-12-29 19:19:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:19:16
	% EndTime: 2019-12-29 19:19:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:19:21
	% EndTime: 2019-12-29 19:19:21
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:19:22
	% EndTime: 2019-12-29 19:19:22
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (44->20), mult. (134->34), div. (0->0), fcn. (92->4), ass. (0->15)
	t134 = sin(qJ(2));
	t136 = cos(qJ(2));
	t146 = r_i_i_C(3) + qJ(3);
	t148 = pkin(2) + r_i_i_C(1);
	t149 = t148 * t134 - t146 * t136;
	t150 = t149 * qJD(2) - t134 * qJD(3);
	t147 = pkin(6) + r_i_i_C(2);
	t135 = sin(qJ(1));
	t145 = qJD(1) * t135;
	t137 = cos(qJ(1));
	t144 = qJD(1) * t137;
	t143 = qJD(2) * t137;
	t141 = -t146 * t134 - t148 * t136;
	t139 = -pkin(1) + t141;
	t1 = [t150 * t135 + (-t147 * t135 + t139 * t137) * qJD(1), (-t146 * t143 + t148 * t145) * t134 + (-t146 * t145 + (-t148 * qJD(2) + qJD(3)) * t137) * t136, -t134 * t145 + t136 * t143, 0, 0; -t150 * t137 + (t139 * t135 + t147 * t137) * qJD(1), -t149 * t144 + (t141 * qJD(2) + qJD(3) * t136) * t135, t135 * qJD(2) * t136 + t134 * t144, 0, 0; 0, -t150, qJD(2) * t134, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:19:21
	% EndTime: 2019-12-29 19:19:21
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (122->38), mult. (376->64), div. (0->0), fcn. (325->6), ass. (0->33)
	t59 = sin(qJ(4));
	t60 = sin(qJ(2));
	t62 = cos(qJ(4));
	t63 = cos(qJ(2));
	t72 = t59 * t63 - t60 * t62;
	t90 = pkin(2) + pkin(3);
	t91 = -qJ(3) * t63 + t90 * t60;
	t94 = t91 * qJD(2) - t60 * qJD(3);
	t71 = t59 * t60 + t62 * t63;
	t93 = qJD(2) - qJD(4);
	t54 = t93 * t71;
	t84 = qJD(2) * t60;
	t92 = t72 * qJD(4) + t62 * t84;
	t61 = sin(qJ(1));
	t86 = qJD(1) * t61;
	t64 = cos(qJ(1));
	t85 = qJD(1) * t64;
	t83 = qJD(2) * t63;
	t82 = qJD(2) * t64;
	t80 = pkin(6) - pkin(7) - r_i_i_C(3);
	t76 = t63 * t82;
	t68 = qJD(1) * t72;
	t49 = t54 * t64 + t61 * t68;
	t67 = qJD(1) * t71;
	t50 = -t59 * t76 + t61 * t67 + t92 * t64;
	t75 = t49 * r_i_i_C(1) + t50 * r_i_i_C(2);
	t51 = -t61 * t54 + t64 * t68;
	t52 = t64 * t67 + (t59 * t83 - t92) * t61;
	t74 = -t51 * r_i_i_C(1) - t52 * r_i_i_C(2);
	t73 = -t93 * t72 * r_i_i_C(1) - t54 * r_i_i_C(2);
	t70 = -qJ(3) * t60 - t90 * t63;
	t66 = -pkin(1) + t70;
	t1 = [-t52 * r_i_i_C(1) + t51 * r_i_i_C(2) + t94 * t61 + (-t80 * t61 + t66 * t64) * qJD(1), (-qJ(3) * t82 + t90 * t86) * t60 + (-qJ(3) * t86 + (-t90 * qJD(2) + qJD(3)) * t64) * t63 - t75, -t60 * t86 + t76, t75, 0; -t50 * r_i_i_C(1) + t49 * r_i_i_C(2) - t94 * t64 + (t66 * t61 + t80 * t64) * qJD(1), -t91 * t85 + (qJD(2) * t70 + qJD(3) * t63) * t61 - t74, t60 * t85 + t61 * t83, t74, 0; 0, -t94 - t73, t84, t73, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:19:19
	% EndTime: 2019-12-29 19:19:20
	% DurationCPUTime: 0.89s
	% Computational Cost: add. (341->73), mult. (1050->123), div. (0->0), fcn. (979->8), ass. (0->52)
	t344 = sin(qJ(4));
	t348 = cos(qJ(2));
	t345 = sin(qJ(2));
	t387 = cos(qJ(4));
	t374 = t345 * t387;
	t396 = -t348 * t344 + t374;
	t390 = pkin(2) + pkin(3);
	t394 = -qJ(3) * t348 + t390 * t345;
	t402 = t394 * qJD(2) - t345 * qJD(3);
	t368 = qJD(4) * t387;
	t369 = qJD(2) * t387;
	t377 = qJD(4) * t344;
	t401 = t348 * t377 + (-t368 + t369) * t345;
	t324 = t345 * t344 + t348 * t387;
	t380 = qJD(2) * t345;
	t318 = t324 * qJD(4) - t344 * t380 - t348 * t369;
	t346 = sin(qJ(1));
	t349 = cos(qJ(1));
	t397 = t396 * qJD(1);
	t316 = t318 * t346 - t397 * t349;
	t379 = qJD(2) * t348;
	t319 = t344 * t379 - t401;
	t354 = qJD(1) * t324;
	t317 = t319 * t346 + t349 * t354;
	t343 = sin(qJ(5));
	t347 = cos(qJ(5));
	t360 = t347 * r_i_i_C(1) - t343 * r_i_i_C(2) + pkin(4);
	t388 = -r_i_i_C(3) - pkin(8);
	t386 = t347 * r_i_i_C(2);
	t395 = qJD(5) * (t343 * r_i_i_C(1) + t386);
	t400 = -t346 * t395 * t396 - t360 * t316 - t388 * t317;
	t378 = qJD(2) * t349;
	t370 = t348 * t378;
	t314 = -t344 * t370 + t346 * t354 + t401 * t349;
	t323 = t324 * t349;
	t383 = t348 * t349;
	t315 = -t349 * t345 * t377 + qJD(2) * t323 - t397 * t346 - t368 * t383;
	t399 = (t344 * t383 - t349 * t374) * t395 + t388 * t314 + t360 * t315;
	t398 = t388 * t318 - t360 * t319 + t324 * t395;
	t389 = pkin(6) - pkin(7);
	t382 = qJD(1) * t346;
	t381 = qJD(1) * t349;
	t376 = qJD(5) * t396;
	t321 = t324 * t346;
	t362 = t321 * t347 + t343 * t349;
	t361 = t321 * t343 - t347 * t349;
	t359 = -qJ(3) * t345 - t390 * t348;
	t355 = -pkin(1) + t359;
	t342 = t343 * t382;
	t313 = -t343 * t381 - t314 * t347 + (-t323 * t343 - t346 * t347) * qJD(5);
	t312 = -t347 * t381 + t314 * t343 + (-t323 * t347 + t343 * t346) * qJD(5);
	t1 = [t342 * r_i_i_C(1) - t360 * t317 + t388 * t316 + (t361 * r_i_i_C(1) + t362 * r_i_i_C(2)) * qJD(5) + t402 * t346 + ((t386 - t389) * t346 + t355 * t349) * qJD(1), (-qJ(3) * t378 + t390 * t382) * t345 + (-qJ(3) * t382 + (-t390 * qJD(2) + qJD(3)) * t349) * t348 - t399, -t345 * t382 + t370, t399, t312 * r_i_i_C(1) - t313 * r_i_i_C(2); -t314 * pkin(4) + t313 * r_i_i_C(1) + t312 * r_i_i_C(2) + t388 * t315 - t402 * t349 + (t355 * t346 + t389 * t349) * qJD(1), -t394 * t381 + (t359 * qJD(2) + qJD(3) * t348) * t346 - t400, t345 * t381 + t346 * t379, t400, (-t317 * t343 - t347 * t382) * r_i_i_C(1) + (-t317 * t347 + t342) * r_i_i_C(2) + (-t362 * r_i_i_C(1) + t361 * r_i_i_C(2)) * qJD(5); 0, -t402 - t398, t380, t398, (t318 * t347 + t343 * t376) * r_i_i_C(2) + (t318 * t343 - t347 * t376) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end