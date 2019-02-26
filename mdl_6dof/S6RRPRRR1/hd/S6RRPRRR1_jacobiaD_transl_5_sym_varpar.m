% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:01
% EndTime: 2019-02-26 21:54:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (276->41), mult. (188->50), div. (0->0), fcn. (122->10), ass. (0->37)
t58 = qJD(2) + qJD(4);
t55 = qJD(5) + t58;
t59 = qJ(2) + pkin(11);
t56 = qJ(4) + t59;
t52 = qJ(5) + t56;
t49 = cos(t52);
t83 = r_i_i_C(2) * t49;
t48 = sin(t52);
t85 = r_i_i_C(1) * t48;
t68 = t83 + t85;
t66 = t68 * t55;
t50 = sin(t56);
t86 = pkin(4) * t50;
t88 = -t58 * t86 - t66;
t70 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t59);
t64 = t70 * qJD(2) + t88;
t81 = t49 * t55;
t75 = r_i_i_C(1) * t81;
t51 = cos(t56);
t80 = t51 * t58;
t87 = -pkin(4) * t80 - t75;
t84 = r_i_i_C(2) * t48;
t82 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3) + pkin(7);
t61 = sin(qJ(1));
t79 = qJD(1) * t61;
t63 = cos(qJ(1));
t78 = qJD(1) * t63;
t74 = t55 * t84;
t72 = qJD(1) * t83;
t73 = t61 * t72 + t63 * t74 + t79 * t85;
t69 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t59);
t71 = t69 * qJD(2) + t87;
t67 = -pkin(4) * t51 - r_i_i_C(1) * t49 - pkin(1) + t69 + t84;
t65 = -t63 * t75 + t73;
t44 = t61 * t74;
t43 = t70 - t86;
t1 = [t63 * qJD(3) - t64 * t61 + (-t82 * t61 + t67 * t63) * qJD(1), -t43 * t79 + t71 * t63 + t73, t78 (t50 * t79 - t63 * t80) * pkin(4) + t65, t65, 0; t61 * qJD(3) + t64 * t63 + (t67 * t61 + t82 * t63) * qJD(1), t44 + t71 * t61 + (t43 - t68) * t78, t79, t44 + t87 * t61 + (-t68 - t86) * t78, -t63 * t72 + t44 + (-t48 * t78 - t61 * t81) * r_i_i_C(1), 0; 0, t64, 0, t88, -t66, 0;];
JaD_transl  = t1;
