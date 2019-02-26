% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function JaD_transl = S6RRRPRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:03
% EndTime: 2019-02-26 22:16:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (149->33), mult. (146->42), div. (0->0), fcn. (95->8), ass. (0->32)
t52 = sin(qJ(2));
t50 = qJD(2) + qJD(3);
t51 = qJ(2) + qJ(3);
t46 = pkin(11) + t51;
t44 = sin(t46);
t45 = cos(t46);
t59 = r_i_i_C(1) * t44 + r_i_i_C(2) * t45;
t47 = sin(t51);
t75 = pkin(3) * t47;
t78 = t59 + t75;
t56 = t78 * t50;
t67 = pkin(2) * qJD(2);
t76 = t52 * t67 + t56;
t48 = cos(t51);
t72 = r_i_i_C(1) * t45;
t77 = -pkin(3) * t48 - t72;
t71 = r_i_i_C(2) * t44;
t69 = r_i_i_C(3) + qJ(4) + pkin(8) + pkin(7);
t68 = t48 * t50;
t53 = sin(qJ(1));
t66 = qJD(1) * t53;
t55 = cos(qJ(1));
t65 = qJD(1) * t55;
t64 = t50 * t72;
t63 = t50 * t71;
t62 = t55 * t63 + t59 * t66;
t54 = cos(qJ(2));
t60 = -pkin(3) * t68 - t54 * t67 - t64;
t58 = -pkin(2) * t54 - pkin(1) + t71 + t77;
t43 = -pkin(2) * t52 - t75;
t38 = t53 * t63;
t1 = [t55 * qJD(4) + t76 * t53 + (-t69 * t53 + t58 * t55) * qJD(1), -t43 * t66 + t60 * t55 + t62, -t55 * t64 + (t47 * t66 - t55 * t68) * pkin(3) + t62, t65, 0, 0; t53 * qJD(4) - t76 * t55 + (t58 * t53 + t69 * t55) * qJD(1), t38 + t60 * t53 + (t43 - t59) * t65, t50 * t53 * t77 - t65 * t78 + t38, t66, 0, 0; 0, -t76, -t56, 0, 0, 0;];
JaD_transl  = t1;
