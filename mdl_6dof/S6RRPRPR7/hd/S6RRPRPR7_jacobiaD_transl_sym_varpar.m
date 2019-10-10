% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR7
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:25
	% EndTime: 2019-10-10 10:15:25
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->20), mult. (134->34), div. (0->0), fcn. (92->4), ass. (0->15)
	t134 = sin(qJ(2));
	t136 = cos(qJ(2));
	t146 = r_i_i_C(3) + qJ(3);
	t148 = pkin(2) + r_i_i_C(1);
	t149 = t148 * t134 - t146 * t136;
	t150 = t149 * qJD(2) - t134 * qJD(3);
	t147 = pkin(7) + r_i_i_C(2);
	t135 = sin(qJ(1));
	t145 = qJD(1) * t135;
	t137 = cos(qJ(1));
	t144 = qJD(1) * t137;
	t143 = qJD(2) * t137;
	t141 = -t146 * t134 - t148 * t136;
	t139 = -pkin(1) + t141;
	t1 = [t150 * t135 + (-t147 * t135 + t139 * t137) * qJD(1), (-t146 * t143 + t148 * t145) * t134 + (-t146 * t145 + (-t148 * qJD(2) + qJD(3)) * t137) * t136, -t134 * t145 + t136 * t143, 0, 0, 0; -t150 * t137 + (t139 * t135 + t147 * t137) * qJD(1), -t149 * t144 + (t141 * qJD(2) + qJD(3) * t136) * t135, t135 * qJD(2) * t136 + t134 * t144, 0, 0, 0; 0, -t150, qJD(2) * t134, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.26s
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
	t92 = qJD(4) * t72 + t62 * t84;
	t61 = sin(qJ(1));
	t86 = qJD(1) * t61;
	t64 = cos(qJ(1));
	t85 = qJD(1) * t64;
	t83 = qJD(2) * t63;
	t82 = qJD(2) * t64;
	t80 = pkin(7) - pkin(8) - r_i_i_C(3);
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
	t1 = [-t52 * r_i_i_C(1) + t51 * r_i_i_C(2) + t94 * t61 + (-t61 * t80 + t64 * t66) * qJD(1), (-qJ(3) * t82 + t90 * t86) * t60 + (-qJ(3) * t86 + (-t90 * qJD(2) + qJD(3)) * t64) * t63 - t75, -t60 * t86 + t76, t75, 0, 0; -t50 * r_i_i_C(1) + t49 * r_i_i_C(2) - t94 * t64 + (t61 * t66 + t64 * t80) * qJD(1), -t91 * t85 + (qJD(2) * t70 + qJD(3) * t63) * t61 - t74, t60 * t85 + t61 * t83, t74, 0, 0; 0, -t94 - t73, t84, t73, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (260->54), mult. (522->76), div. (0->0), fcn. (436->8), ass. (0->43)
	t69 = qJ(4) + pkin(10);
	t67 = sin(t69);
	t68 = cos(t69);
	t72 = sin(qJ(2));
	t75 = cos(qJ(2));
	t90 = t67 * t75 - t68 * t72;
	t74 = cos(qJ(4));
	t108 = pkin(4) * t74 + pkin(2) + pkin(3);
	t71 = sin(qJ(4));
	t94 = pkin(4) * t71 + qJ(3);
	t82 = t108 * t72 - t94 * t75;
	t112 = -t82 * qJD(2) + t72 * qJD(3);
	t110 = qJD(2) - qJD(4);
	t89 = t67 * t72 + t68 * t75;
	t61 = t110 * t89;
	t111 = (pkin(7) - r_i_i_C(3) - qJ(5) - pkin(8)) * qJD(1);
	t102 = qJD(2) * t72;
	t109 = t90 * qJD(4) + t68 * t102;
	t105 = pkin(4) * qJD(4);
	t73 = sin(qJ(1));
	t104 = qJD(1) * t73;
	t76 = cos(qJ(1));
	t103 = qJD(1) * t76;
	t101 = qJD(2) * t75;
	t97 = t76 * t101;
	t86 = qJD(1) * t90;
	t56 = t61 * t76 + t73 * t86;
	t85 = qJD(1) * t89;
	t57 = t109 * t76 - t67 * t97 + t73 * t85;
	t93 = t56 * r_i_i_C(1) + t57 * r_i_i_C(2);
	t58 = -t73 * t61 + t76 * t86;
	t59 = t76 * t85 + (t67 * t101 - t109) * t73;
	t92 = -t58 * r_i_i_C(1) - t59 * r_i_i_C(2);
	t91 = -t110 * t90 * r_i_i_C(1) - t61 * r_i_i_C(2);
	t88 = t71 * t75 - t72 * t74;
	t87 = t71 * t72 + t74 * t75;
	t84 = t88 * qJD(4);
	t83 = -t108 * t75 - t94 * t72;
	t80 = -qJD(5) + (-pkin(1) + t83) * qJD(1);
	t79 = t110 * t87;
	t78 = t83 * qJD(2) + qJD(3) * t75 + t87 * t105;
	t77 = -t88 * t105 + t112;
	t1 = [-t59 * r_i_i_C(1) + t58 * r_i_i_C(2) + t80 * t76 + (pkin(4) * t84 - t111 - t112) * t73, t82 * t104 + t78 * t76 - t93, -t72 * t104 + t97, (t88 * t104 + t79 * t76) * pkin(4) + t93, -t103, 0; -t57 * r_i_i_C(1) + t56 * r_i_i_C(2) + t80 * t73 + (t77 + t111) * t76, -t103 * t82 + t78 * t73 - t92, t73 * t101 + t72 * t103, (-t88 * t103 + t79 * t73) * pkin(4) + t92, -t104, 0; 0, t77 - t91, t102, (-t88 * qJD(2) + t84) * pkin(4) + t91, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:26
	% EndTime: 2019-10-10 10:15:27
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (687->90), mult. (1196->138), div. (0->0), fcn. (1090->10), ass. (0->63)
	t398 = qJ(4) + pkin(10);
	t355 = sin(t398);
	t363 = cos(qJ(2));
	t359 = sin(qJ(2));
	t392 = cos(t398);
	t388 = t359 * t392;
	t336 = -t363 * t355 + t388;
	t358 = sin(qJ(4));
	t393 = pkin(4) * t358 + qJ(3);
	t362 = cos(qJ(4));
	t412 = t362 * pkin(4) + pkin(2) + pkin(3);
	t374 = t412 * t359 - t393 * t363;
	t424 = t374 * qJD(2) - t359 * qJD(3);
	t386 = qJD(4) * t392;
	t387 = qJD(2) * t392;
	t402 = qJD(4) * t355;
	t423 = t363 * t402 + (-t386 + t387) * t359;
	t335 = t359 * t355 + t363 * t392;
	t404 = qJD(2) * t359;
	t329 = t335 * qJD(4) - t355 * t404 - t363 * t387;
	t360 = sin(qJ(1));
	t364 = cos(qJ(1));
	t419 = t336 * qJD(1);
	t327 = t329 * t360 - t419 * t364;
	t403 = qJD(2) * t363;
	t330 = t355 * t403 - t423;
	t372 = qJD(1) * t335;
	t328 = t330 * t360 + t364 * t372;
	t357 = sin(qJ(6));
	t361 = cos(qJ(6));
	t390 = -t361 * r_i_i_C(1) + t357 * r_i_i_C(2);
	t383 = pkin(5) - t390;
	t414 = -r_i_i_C(3) - pkin(9);
	t413 = t361 * r_i_i_C(2);
	t418 = qJD(6) * (t357 * r_i_i_C(1) + t413);
	t422 = -t336 * t360 * t418 - t383 * t327 - t414 * t328;
	t394 = t364 * t403;
	t325 = -t355 * t394 + t360 * t372 + t364 * t423;
	t334 = t335 * t364;
	t407 = t363 * t364;
	t326 = -t364 * t359 * t402 + qJD(2) * t334 - t419 * t360 - t386 * t407;
	t421 = (t355 * t407 - t364 * t388) * t418 + t414 * t325 + t383 * t326;
	t420 = t414 * t329 - t383 * t330 + t335 * t418;
	t411 = pkin(7) - qJ(5) - pkin(8);
	t410 = pkin(4) * qJD(4);
	t409 = t328 * t357;
	t406 = qJD(1) * t360;
	t405 = qJD(1) * t364;
	t401 = qJD(6) * t357;
	t400 = qJD(6) * t361;
	t391 = -t328 * t361 + t357 * t406;
	t385 = t358 * t363 - t359 * t362;
	t384 = t358 * t359 + t362 * t363;
	t376 = t385 * qJD(4);
	t375 = -t393 * t359 - t412 * t363;
	t368 = -qJD(5) + (-pkin(1) + t375) * qJD(1);
	t367 = (qJD(2) - qJD(4)) * t384;
	t366 = t375 * qJD(2) + qJD(3) * t363 + t384 * t410;
	t365 = -t385 * t410 - t424;
	t332 = t335 * t360;
	t324 = -t357 * t405 - t325 * t361 + (-t334 * t357 - t360 * t361) * qJD(6);
	t323 = -t361 * t405 + t325 * t357 + (-t334 * t361 + t357 * t360) * qJD(6);
	t1 = [(t332 * t401 + t391) * r_i_i_C(1) + (t332 * t400 + t409) * r_i_i_C(2) - t328 * pkin(5) + t414 * t327 + (t390 * qJD(6) + t368) * t364 + (pkin(4) * t376 + (-t411 + t413) * qJD(1) + t424) * t360, t366 * t364 + t374 * t406 - t421, -t359 * t406 + t394, (t367 * t364 + t385 * t406) * pkin(4) + t421, -t405, t323 * r_i_i_C(1) - t324 * r_i_i_C(2); -t325 * pkin(5) + t324 * r_i_i_C(1) + t323 * r_i_i_C(2) + t414 * t326 + t368 * t360 + (t411 * qJD(1) + t365) * t364, t366 * t360 - t374 * t405 - t422, t359 * t405 + t360 * t403, (t367 * t360 - t385 * t405) * pkin(4) + t422, -t406, (-t361 * t406 - t409) * r_i_i_C(1) + t391 * r_i_i_C(2) + ((-t332 * t361 - t357 * t364) * r_i_i_C(1) + (t332 * t357 - t361 * t364) * r_i_i_C(2)) * qJD(6); 0, t365 - t420, t404, (-t385 * qJD(2) + t376) * pkin(4) + t420, 0, (t329 * t361 + t336 * t401) * r_i_i_C(2) + (t329 * t357 - t336 * t400) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end