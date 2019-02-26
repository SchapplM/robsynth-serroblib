% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:03
% EndTime: 2019-02-26 22:16:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (294->41), mult. (196->49), div. (0->0), fcn. (127->10), ass. (0->35)
t63 = sin(qJ(2));
t62 = qJ(2) + qJ(3);
t57 = pkin(11) + t62;
t48 = -pkin(3) * sin(t62) - pkin(4) * sin(t57);
t61 = qJD(2) + qJD(3);
t56 = qJD(5) + t61;
t55 = qJ(5) + t57;
t52 = cos(t55);
t85 = r_i_i_C(2) * t52;
t51 = sin(t55);
t87 = r_i_i_C(1) * t51;
t71 = t85 + t87;
t68 = t71 * t56;
t67 = t48 * t61 - t68;
t82 = pkin(2) * qJD(2);
t88 = -t63 * t82 + t67;
t72 = -pkin(3) * cos(t62) - pkin(4) * cos(t57);
t83 = t52 * t56;
t79 = r_i_i_C(1) * t83;
t73 = t72 * t61 - t79;
t86 = r_i_i_C(2) * t51;
t84 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8) + pkin(7);
t64 = sin(qJ(1));
t81 = qJD(1) * t64;
t66 = cos(qJ(1));
t80 = qJD(1) * t66;
t78 = t56 * t86;
t76 = qJD(1) * t85;
t77 = t64 * t76 + t66 * t78 + t81 * t87;
t65 = cos(qJ(2));
t74 = -t65 * t82 + t73;
t70 = -t65 * pkin(2) - r_i_i_C(1) * t52 - pkin(1) + t72 + t86;
t46 = t64 * t78;
t45 = -t63 * pkin(2) + t48;
t1 = [t66 * qJD(4) - t88 * t64 + (-t84 * t64 + t70 * t66) * qJD(1), -t45 * t81 + t74 * t66 + t77, -t48 * t81 + t73 * t66 + t77, t80, -t66 * t79 + t77, 0; t64 * qJD(4) + t88 * t66 + (t70 * t64 + t84 * t66) * qJD(1), t46 + t74 * t64 + (t45 - t71) * t80, t46 + t73 * t64 + (t48 - t71) * t80, t81, -t66 * t76 + t46 + (-t51 * t80 - t64 * t83) * r_i_i_C(1), 0; 0, t88, t67, 0, -t68, 0;];
JaD_transl  = t1;
