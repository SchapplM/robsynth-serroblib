% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:43
% EndTime: 2019-02-26 20:34:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (166->33), mult. (124->42), div. (0->0), fcn. (81->9), ass. (0->29)
t50 = qJD(4) + qJD(5);
t49 = pkin(11) + qJ(4);
t47 = qJ(5) + t49;
t42 = cos(t47);
t66 = r_i_i_C(2) * t42;
t41 = sin(t47);
t68 = r_i_i_C(1) * t41;
t56 = t66 + t68;
t54 = t56 * t50;
t43 = sin(t49);
t69 = pkin(4) * t43;
t70 = qJD(4) * t69 + t54;
t67 = r_i_i_C(2) * t41;
t65 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(3);
t64 = t42 * t50;
t51 = qJ(1) + pkin(10);
t44 = sin(t51);
t63 = qJD(1) * t44;
t46 = cos(t51);
t62 = qJD(1) * t46;
t45 = cos(t49);
t61 = qJD(4) * t45;
t60 = r_i_i_C(1) * t64;
t59 = t50 * t67;
t57 = qJD(1) * t66;
t55 = -r_i_i_C(1) * t42 - pkin(4) * t45 - cos(pkin(11)) * pkin(3) - pkin(2) + t67;
t53 = t44 * t57 + t63 * t68 + (t59 - t60) * t46;
t36 = t44 * t59;
t1 = [t46 * qJD(3) + t70 * t44 + (-cos(qJ(1)) * pkin(1) - t65 * t44 + t55 * t46) * qJD(1), 0, t62 (t43 * t63 - t46 * t61) * pkin(4) + t53, t53, 0; t44 * qJD(3) - t70 * t46 + (-sin(qJ(1)) * pkin(1) + t65 * t46 + t55 * t44) * qJD(1), 0, t63, t36 + (-pkin(4) * t61 - t60) * t44 + (-t56 - t69) * t62, -t46 * t57 + t36 + (-t41 * t62 - t44 * t64) * r_i_i_C(1), 0; 0, 0, 0, -t70, -t54, 0;];
JaD_transl  = t1;
