% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:25
% DurationCPUTime: 0.20s
% Computational Cost: add. (234->45), mult. (424->68), div. (0->0), fcn. (364->8), ass. (0->35)
t67 = pkin(10) + qJ(5);
t65 = sin(t67);
t66 = cos(t67);
t70 = sin(qJ(2));
t72 = cos(qJ(2));
t81 = t65 * t72 - t66 * t70;
t85 = pkin(4) * sin(pkin(10)) + qJ(3);
t99 = pkin(2) + cos(pkin(10)) * pkin(4) + pkin(3);
t100 = t70 * t99 - t72 * t85;
t103 = qJD(2) * t100 - t70 * qJD(3);
t102 = qJD(2) - qJD(5);
t80 = t65 * t70 + t66 * t72;
t59 = t102 * t80;
t94 = qJD(2) * t70;
t101 = qJD(5) * t81 + t66 * t94;
t71 = sin(qJ(1));
t96 = qJD(1) * t71;
t73 = cos(qJ(1));
t95 = qJD(1) * t73;
t93 = qJD(2) * t72;
t92 = qJD(2) * t73;
t90 = pkin(7) - r_i_i_C(3) - pkin(8) - qJ(4);
t86 = t72 * t92;
t79 = qJD(1) * t81;
t54 = t59 * t73 + t71 * t79;
t78 = qJD(1) * t80;
t55 = t101 * t73 - t65 * t86 + t71 * t78;
t84 = t54 * r_i_i_C(1) + t55 * r_i_i_C(2);
t56 = -t59 * t71 + t73 * t79;
t57 = t73 * t78 + (t65 * t93 - t101) * t71;
t83 = -t56 * r_i_i_C(1) - t57 * r_i_i_C(2);
t82 = -r_i_i_C(1) * t102 * t81 - t59 * r_i_i_C(2);
t77 = -t70 * t85 - t72 * t99;
t75 = -pkin(1) + t77;
t1 = [-t57 * r_i_i_C(1) + t56 * r_i_i_C(2) - t73 * qJD(4) + t103 * t71 + (-t71 * t90 + t73 * t75) * qJD(1) (-t85 * t92 + t96 * t99) * t70 + (-t85 * t96 + (-qJD(2) * t99 + qJD(3)) * t73) * t72 - t84, -t70 * t96 + t86, -t95, t84, 0; -t55 * r_i_i_C(1) + t54 * r_i_i_C(2) - t71 * qJD(4) - t103 * t73 + (t71 * t75 + t73 * t90) * qJD(1), -t100 * t95 + (qJD(2) * t77 + qJD(3) * t72) * t71 - t83, t70 * t95 + t71 * t93, -t96, t83, 0; 0, -t103 - t82, t94, 0, t82, 0;];
JaD_transl  = t1;
