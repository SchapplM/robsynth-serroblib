% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:32
% EndTime: 2019-02-26 22:04:32
% DurationCPUTime: 0.14s
% Computational Cost: add. (221->46), mult. (265->57), div. (0->0), fcn. (180->6), ass. (0->36)
t68 = qJ(2) + qJ(3);
t66 = cos(t68);
t97 = r_i_i_C(1) + qJ(4);
t109 = t97 * t66;
t65 = sin(t68);
t63 = t65 * qJD(4);
t67 = qJD(2) + qJD(3);
t69 = sin(qJ(2));
t96 = pkin(2) * qJD(2);
t88 = t69 * t96;
t110 = pkin(3) + pkin(4);
t92 = -r_i_i_C(2) + t110;
t111 = (-t92 * t65 + t109) * t67 - (r_i_i_C(3) + qJ(5) - pkin(8) - pkin(7)) * qJD(1) + t63 - t88;
t70 = sin(qJ(1));
t100 = t66 * t70;
t72 = cos(qJ(1));
t94 = qJD(1) * t72;
t108 = t67 * t100 + t65 * t94;
t102 = pkin(2) * t69;
t99 = t67 * t65;
t98 = t67 * t72;
t95 = qJD(1) * t70;
t93 = qJD(4) * t66;
t90 = t66 * t98;
t87 = t110 * t67;
t86 = t110 * t72;
t84 = t65 * t95;
t82 = t97 * t65;
t80 = t97 * t70;
t79 = r_i_i_C(2) * t90 + t110 * t84 + t72 * t93;
t78 = t108 * r_i_i_C(2) + t109 * t94 + t70 * t93;
t76 = r_i_i_C(2) * t99 + t109 * t67 - t65 * t87 + t63;
t71 = cos(qJ(2));
t75 = -qJD(5) + (-t71 * pkin(2) - t92 * t66 - pkin(1) - t82) * qJD(1);
t74 = -t71 * t96 + (-t110 * t66 - t82) * t67;
t1 = [-t111 * t70 + t75 * t72 (-r_i_i_C(2) * t65 + t102 - t109) * t95 + t74 * t72 + t79 (-r_i_i_C(2) * t95 - t97 * t98) * t65 + (-qJD(1) * t80 - t67 * t86) * t66 + t79, -t84 + t90, -t94, 0; t111 * t72 + t75 * t70 (-t110 * t65 - t102) * t94 + t74 * t70 + t78, -t87 * t100 + (-qJD(1) * t86 - t67 * t80) * t65 + t78, t108, -t95, 0; 0, t76 - t88, t76, t99, 0, 0;];
JaD_transl  = t1;
