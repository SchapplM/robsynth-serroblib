% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:37:45
% EndTime: 2019-02-26 19:37:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (196->37), mult. (166->50), div. (0->0), fcn. (107->8), ass. (0->37)
t57 = qJ(2) + qJ(3);
t55 = qJ(4) + t57;
t51 = cos(t55);
t56 = qJD(2) + qJD(3);
t52 = qJD(4) + t56;
t76 = t51 * t52;
t69 = r_i_i_C(1) * t76;
t54 = cos(t57);
t75 = t54 * t56;
t82 = -pkin(3) * t75 - t69;
t53 = sin(t57);
t80 = pkin(3) * t53;
t49 = t56 * t80;
t58 = sin(qJ(2));
t73 = pkin(2) * qJD(2);
t39 = -t58 * t73 - t49;
t78 = r_i_i_C(2) * t51;
t50 = sin(t55);
t79 = r_i_i_C(1) * t50;
t64 = t78 + t79;
t81 = t64 * t52 - t39;
t77 = t50 * t52;
t74 = r_i_i_C(1) * t77 + r_i_i_C(2) * t76;
t59 = sin(qJ(1));
t72 = qJD(1) * t59;
t61 = cos(qJ(1));
t71 = qJD(1) * t61;
t68 = r_i_i_C(2) * t77;
t66 = qJD(1) * t78;
t67 = t59 * t66 + t61 * t68 + t72 * t79;
t60 = cos(qJ(2));
t65 = -t60 * t73 + t82;
t63 = -pkin(2) * t60 - pkin(3) * t54 - r_i_i_C(1) * t51 + r_i_i_C(2) * t50 - pkin(1);
t62 = -t61 * t69 + t67;
t48 = -pkin(2) * t58 - t80;
t41 = t59 * t68;
t1 = [t81 * t59 + (r_i_i_C(3) * t59 + t63 * t61) * qJD(1), -t48 * t72 + t65 * t61 + t67 (t53 * t72 - t61 * t75) * pkin(3) + t62, t62, 0; -t81 * t61 + (-r_i_i_C(3) * t61 + t63 * t59) * qJD(1), t41 + t65 * t59 + (t48 - t64) * t71, t41 + t82 * t59 + (-t64 - t80) * t71, -t61 * t66 + t41 + (-t50 * t71 - t59 * t76) * r_i_i_C(1), 0; 0, -t39 + t74, t49 + t74, t74, 0;];
JaD_transl  = t1;
