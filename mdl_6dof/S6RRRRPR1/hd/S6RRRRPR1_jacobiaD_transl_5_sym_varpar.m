% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:46
% EndTime: 2019-02-26 22:30:46
% DurationCPUTime: 0.21s
% Computational Cost: add. (323->43), mult. (210->51), div. (0->0), fcn. (136->10), ass. (0->40)
t64 = sin(qJ(2));
t63 = qJ(2) + qJ(3);
t58 = sin(t63);
t62 = qJD(2) + qJD(3);
t57 = qJD(4) + t62;
t61 = qJ(4) + t63;
t54 = pkin(11) + t61;
t52 = sin(t54);
t53 = cos(t54);
t72 = r_i_i_C(1) * t52 + r_i_i_C(2) * t53;
t55 = sin(t61);
t90 = pkin(4) * t55;
t94 = t72 + t90;
t69 = t94 * t57;
t91 = pkin(3) * t62;
t68 = -t58 * t91 - t69;
t82 = pkin(2) * qJD(2);
t95 = -t64 * t82 + t68;
t59 = cos(t63);
t87 = r_i_i_C(1) * t53;
t78 = t57 * t87;
t56 = cos(t61);
t83 = t56 * t57;
t73 = -pkin(4) * t83 - t59 * t91 - t78;
t93 = -pkin(4) * t56 - t87;
t86 = r_i_i_C(2) * t52;
t84 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8) + pkin(7);
t65 = sin(qJ(1));
t81 = qJD(1) * t65;
t67 = cos(qJ(1));
t80 = qJD(1) * t67;
t77 = t57 * t86;
t76 = t67 * t77 + t72 * t81;
t66 = cos(qJ(2));
t74 = -t66 * t82 + t73;
t49 = -pkin(3) * t58 - t90;
t71 = -t66 * pkin(2) - pkin(3) * t59 - pkin(1) + t86 + t93;
t47 = t65 * t77;
t46 = -t64 * pkin(2) + t49;
t1 = [t67 * qJD(5) - t95 * t65 + (-t84 * t65 + t71 * t67) * qJD(1), -t46 * t81 + t74 * t67 + t76, -t49 * t81 + t73 * t67 + t76, -t67 * t78 + (t55 * t81 - t67 * t83) * pkin(4) + t76, t80, 0; t65 * qJD(5) + t95 * t67 + (t71 * t65 + t84 * t67) * qJD(1), t47 + t74 * t65 + (t46 - t72) * t80, t47 + t73 * t65 + (t49 - t72) * t80, t93 * t65 * t57 - t80 * t94 + t47, t81, 0; 0, t95, t68, -t69, 0, 0;];
JaD_transl  = t1;
