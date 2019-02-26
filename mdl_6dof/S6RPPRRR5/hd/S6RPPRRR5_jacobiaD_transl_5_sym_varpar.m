% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:09
% EndTime: 2019-02-26 20:37:09
% DurationCPUTime: 0.10s
% Computational Cost: add. (90->30), mult. (134->39), div. (0->0), fcn. (87->6), ass. (0->28)
t43 = cos(qJ(4));
t39 = qJD(4) + qJD(5);
t40 = qJ(4) + qJ(5);
t37 = sin(t40);
t59 = r_i_i_C(2) * t37;
t38 = cos(t40);
t60 = r_i_i_C(1) * t38;
t51 = (t59 - t60) * t39;
t55 = pkin(4) * qJD(4);
t65 = -t43 * t55 + t51;
t63 = qJD(3) - t65;
t62 = pkin(4) * t43;
t61 = r_i_i_C(1) * t37;
t58 = r_i_i_C(2) * t38;
t42 = sin(qJ(1));
t57 = t39 * t42;
t44 = cos(qJ(1));
t56 = t39 * t44;
t54 = qJD(1) * t42;
t36 = qJD(1) * t44;
t53 = r_i_i_C(3) - qJ(2) + pkin(8) + pkin(7);
t49 = -t58 - t61;
t41 = sin(qJ(4));
t47 = -pkin(4) * t41 - pkin(1) - qJ(3) + t49;
t46 = t49 * t39 - t41 * t55;
t34 = t36 * t60;
t33 = t54 * t59;
t1 = [t44 * qJD(2) - t63 * t42 + (t53 * t42 + t47 * t44) * qJD(1), t36, -t54, t33 + (-t60 - t62) * t54 + t46 * t44, -t56 * t58 + t33 + (-t37 * t56 - t38 * t54) * r_i_i_C(1), 0; t42 * qJD(2) + t63 * t44 + (t47 * t42 - t53 * t44) * qJD(1), t54, t36, t34 + (-t59 + t62) * t36 + t46 * t42, -t57 * t61 + t34 + (-t37 * t36 - t38 * t57) * r_i_i_C(2), 0; 0, 0, 0, t65, t51, 0;];
JaD_transl  = t1;
