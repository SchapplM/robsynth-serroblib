% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:50
% EndTime: 2019-02-26 21:00:50
% DurationCPUTime: 0.17s
% Computational Cost: add. (195->36), mult. (150->44), div. (0->0), fcn. (97->10), ass. (0->33)
t59 = sin(qJ(3));
t56 = qJD(3) + qJD(4);
t58 = qJ(3) + qJ(4);
t52 = pkin(11) + t58;
t48 = sin(t52);
t49 = cos(t52);
t64 = r_i_i_C(1) * t48 + r_i_i_C(2) * t49;
t53 = sin(t58);
t80 = pkin(4) * t53;
t83 = t64 + t80;
t61 = t83 * t56;
t72 = pkin(3) * qJD(3);
t81 = t59 * t72 + t61;
t54 = cos(t58);
t77 = r_i_i_C(1) * t49;
t82 = -pkin(4) * t54 - t77;
t76 = r_i_i_C(2) * t48;
t74 = r_i_i_C(3) + qJ(5) + pkin(8) + pkin(7);
t73 = t54 * t56;
t57 = qJ(1) + pkin(10);
t50 = sin(t57);
t71 = qJD(1) * t50;
t51 = cos(t57);
t70 = qJD(1) * t51;
t69 = t56 * t77;
t68 = t56 * t76;
t67 = t51 * t68 + t64 * t71;
t60 = cos(qJ(3));
t65 = -pkin(4) * t73 - t60 * t72 - t69;
t63 = -pkin(3) * t60 - pkin(2) + t76 + t82;
t47 = -pkin(3) * t59 - t80;
t40 = t50 * t68;
t1 = [t51 * qJD(5) + t81 * t50 + (-cos(qJ(1)) * pkin(1) - t74 * t50 + t63 * t51) * qJD(1), 0, -t47 * t71 + t65 * t51 + t67, -t51 * t69 + (-t51 * t73 + t53 * t71) * pkin(4) + t67, t70, 0; t50 * qJD(5) - t81 * t51 + (-sin(qJ(1)) * pkin(1) + t74 * t51 + t63 * t50) * qJD(1), 0, t40 + t65 * t50 + (t47 - t64) * t70, t50 * t56 * t82 - t70 * t83 + t40, t71, 0; 0, 0, -t81, -t61, 0, 0;];
JaD_transl  = t1;
