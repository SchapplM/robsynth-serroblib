% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR8_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR8_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:40
% EndTime: 2019-02-26 21:18:40
% DurationCPUTime: 0.08s
% Computational Cost: add. (85->27), mult. (126->41), div. (0->0), fcn. (81->6), ass. (0->27)
t43 = cos(qJ(3));
t62 = pkin(3) * t43;
t40 = qJ(3) + qJ(4);
t38 = cos(t40);
t61 = r_i_i_C(1) * t38;
t37 = sin(t40);
t60 = r_i_i_C(2) * t37;
t39 = qJD(3) + qJD(4);
t59 = t37 * t39;
t58 = t38 * t39;
t42 = sin(qJ(1));
t57 = qJD(1) * t42;
t44 = cos(qJ(1));
t56 = qJD(1) * t44;
t41 = sin(qJ(3));
t55 = qJD(3) * t41;
t54 = -pkin(1) - r_i_i_C(3) - pkin(8) - pkin(7);
t53 = r_i_i_C(1) * t59;
t52 = qJD(1) * t61;
t51 = qJD(3) * t62;
t50 = -r_i_i_C(1) * t58 + r_i_i_C(2) * t59;
t49 = -r_i_i_C(1) * t37 - r_i_i_C(2) * t38;
t48 = t42 * t52 - t57 * t60 + (t58 * r_i_i_C(2) + t53) * t44;
t47 = pkin(3) * t41 + qJ(2) - t49;
t46 = t51 + qJD(2) + (-t60 + t61) * t39;
t35 = t44 * t52;
t1 = [t46 * t44 + (-t47 * t42 + t54 * t44) * qJD(1), t56, t35 + (-t60 + t62) * t56 + (-pkin(3) * t55 + t49 * t39) * t42, -t42 * t53 + t35 + (-t37 * t56 - t42 * t58) * r_i_i_C(2), 0, 0; t46 * t42 + (t54 * t42 + t47 * t44) * qJD(1), t57 (t43 * t57 + t44 * t55) * pkin(3) + t48, t48, 0, 0; 0, 0, t50 - t51, t50, 0, 0;];
JaD_transl  = t1;
