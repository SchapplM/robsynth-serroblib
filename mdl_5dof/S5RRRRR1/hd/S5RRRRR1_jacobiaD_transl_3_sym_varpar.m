% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
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

function JaD_transl = S5RRRRR1_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR1_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobiaD_transl_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:37:45
% EndTime: 2019-02-26 19:37:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (77->25), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
t41 = qJD(2) + qJD(3);
t42 = qJ(2) + qJ(3);
t40 = cos(t42);
t60 = r_i_i_C(2) * t40;
t39 = sin(t42);
t61 = r_i_i_C(1) * t39;
t49 = t60 + t61;
t43 = sin(qJ(2));
t62 = pkin(2) * t43;
t50 = qJD(2) * t62;
t63 = t49 * t41 + t50;
t59 = t39 * t41;
t58 = t40 * t41;
t57 = r_i_i_C(1) * t59 + r_i_i_C(2) * t58;
t44 = sin(qJ(1));
t56 = qJD(1) * t44;
t46 = cos(qJ(1));
t55 = qJD(1) * t46;
t45 = cos(qJ(2));
t54 = qJD(2) * t45;
t53 = r_i_i_C(1) * t58;
t52 = r_i_i_C(2) * t59;
t51 = qJD(1) * t60;
t48 = -t45 * pkin(2) - r_i_i_C(1) * t40 + r_i_i_C(2) * t39 - pkin(1);
t47 = t44 * t51 + t56 * t61 + (t52 - t53) * t46;
t32 = t44 * t52;
t1 = [t63 * t44 + (r_i_i_C(3) * t44 + t48 * t46) * qJD(1) (t43 * t56 - t46 * t54) * pkin(2) + t47, t47, 0, 0; -t63 * t46 + (-r_i_i_C(3) * t46 + t48 * t44) * qJD(1), t32 + (-pkin(2) * t54 - t53) * t44 + (-t49 - t62) * t55, -t46 * t51 + t32 + (-t39 * t55 - t44 * t58) * r_i_i_C(1), 0, 0; 0, t50 + t57, t57, 0, 0;];
JaD_transl  = t1;
