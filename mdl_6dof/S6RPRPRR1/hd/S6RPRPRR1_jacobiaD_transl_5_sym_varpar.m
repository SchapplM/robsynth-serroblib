% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:03
% EndTime: 2019-02-26 20:49:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (173->33), mult. (136->42), div. (0->0), fcn. (88->10), ass. (0->28)
t55 = qJ(3) + pkin(11);
t45 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t55);
t54 = qJD(3) + qJD(5);
t52 = qJ(5) + t55;
t47 = cos(t52);
t72 = r_i_i_C(2) * t47;
t46 = sin(t52);
t74 = r_i_i_C(1) * t46;
t61 = t72 + t74;
t59 = t61 * t54;
t75 = t45 * qJD(3) - t59;
t73 = r_i_i_C(2) * t46;
t71 = r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
t70 = t47 * t54;
t56 = qJ(1) + pkin(10);
t49 = sin(t56);
t69 = qJD(1) * t49;
t51 = cos(t56);
t68 = qJD(1) * t51;
t67 = r_i_i_C(1) * t70;
t66 = t54 * t73;
t64 = qJD(1) * t72;
t65 = t49 * t64 + t51 * t66 + t69 * t74;
t62 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t55);
t63 = t62 * qJD(3) - t67;
t60 = -r_i_i_C(1) * t47 - pkin(2) + t62 + t73;
t38 = t49 * t66;
t1 = [t51 * qJD(4) - t75 * t49 + (-cos(qJ(1)) * pkin(1) - t71 * t49 + t60 * t51) * qJD(1), 0, -t45 * t69 + t63 * t51 + t65, t68, -t51 * t67 + t65, 0; t49 * qJD(4) + t75 * t51 + (-sin(qJ(1)) * pkin(1) + t71 * t51 + t60 * t49) * qJD(1), 0, t38 + t63 * t49 + (t45 - t61) * t68, t69, -t51 * t64 + t38 + (-t46 * t68 - t49 * t70) * r_i_i_C(1), 0; 0, 0, t75, 0, -t59, 0;];
JaD_transl  = t1;
