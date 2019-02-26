% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:51
% EndTime: 2019-02-26 21:14:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (252->39), mult. (174->49), div. (0->0), fcn. (111->10), ass. (0->38)
t61 = sin(qJ(3));
t57 = qJD(3) + qJD(4);
t53 = qJD(5) + t57;
t60 = qJ(3) + qJ(4);
t56 = qJ(5) + t60;
t50 = cos(t56);
t82 = r_i_i_C(2) * t50;
t49 = sin(t56);
t84 = r_i_i_C(1) * t49;
t67 = t82 + t84;
t65 = t67 * t53;
t54 = sin(t60);
t85 = pkin(4) * t54;
t63 = -t57 * t85 - t65;
t78 = pkin(3) * qJD(3);
t87 = -t61 * t78 + t63;
t80 = t50 * t53;
t73 = r_i_i_C(1) * t80;
t55 = cos(t60);
t79 = t55 * t57;
t86 = -pkin(4) * t79 - t73;
t83 = r_i_i_C(2) * t49;
t81 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
t58 = qJ(1) + pkin(11);
t51 = sin(t58);
t77 = qJD(1) * t51;
t52 = cos(t58);
t76 = qJD(1) * t52;
t72 = t53 * t83;
t70 = qJD(1) * t82;
t71 = t51 * t70 + t52 * t72 + t77 * t84;
t62 = cos(qJ(3));
t68 = -t62 * t78 + t86;
t66 = -t62 * pkin(3) - pkin(4) * t55 - r_i_i_C(1) * t50 - pkin(2) + t83;
t64 = -t52 * t73 + t71;
t48 = -t61 * pkin(3) - t85;
t41 = t51 * t72;
t1 = [-t87 * t51 + (-cos(qJ(1)) * pkin(1) - t81 * t51 + t66 * t52) * qJD(1), 0, -t48 * t77 + t68 * t52 + t71 (-t52 * t79 + t54 * t77) * pkin(4) + t64, t64, 0; t87 * t52 + (-sin(qJ(1)) * pkin(1) + t81 * t52 + t66 * t51) * qJD(1), 0, t41 + t68 * t51 + (t48 - t67) * t76, t41 + t86 * t51 + (-t67 - t85) * t76, -t52 * t70 + t41 + (-t49 * t76 - t51 * t80) * r_i_i_C(1), 0; 0, 0, t87, t63, -t65, 0;];
JaD_transl  = t1;
