% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR6_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR6_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (77->42), mult. (230->67), div. (0->0), fcn. (185->6), ass. (0->27)
t36 = sin(qJ(2));
t50 = t36 * qJD(3);
t38 = cos(qJ(2));
t55 = pkin(2) + pkin(3);
t56 = -t38 * qJ(3) + t55 * t36;
t57 = t56 * qJD(2) - t50;
t37 = sin(qJ(1));
t53 = qJD(1) * t37;
t39 = cos(qJ(1));
t52 = qJD(1) * t39;
t51 = qJD(2) * t39;
t49 = pkin(7) - r_i_i_C(3) - qJ(4);
t34 = sin(pkin(10));
t35 = cos(pkin(10));
t48 = t34 * t38 - t35 * t36;
t47 = t34 * t36 + t35 * t38;
t46 = -qJ(3) * t36 - t55 * t38;
t44 = t47 * t37;
t43 = t47 * t39;
t42 = qJD(1) * t48;
t41 = t48 * qJD(2);
t40 = -pkin(1) + t46;
t33 = qJD(1) * t43 + t37 * t41;
t32 = -qJD(2) * t44 + t39 * t42;
t31 = -qJD(1) * t44 + t39 * t41;
t30 = qJD(2) * t43 + t37 * t42;
t1 = [-t33 * r_i_i_C(1) + t32 * r_i_i_C(2) - t39 * qJD(4) + t57 * t37 + (-t49 * t37 + t40 * t39) * qJD(1), -t30 * r_i_i_C(1) + t31 * r_i_i_C(2) + (-qJ(3) * t51 + t55 * t53) * t36 + (-qJ(3) * t53 + (-t55 * qJD(2) + qJD(3)) * t39) * t38, -t36 * t53 + t38 * t51, -t52, 0, 0; t31 * r_i_i_C(1) + t30 * r_i_i_C(2) - t37 * qJD(4) - t57 * t39 + (t40 * t37 + t49 * t39) * qJD(1), t32 * r_i_i_C(1) + t33 * r_i_i_C(2) - t56 * t52 + (t46 * qJD(2) + qJD(3) * t38) * t37, t37 * qJD(2) * t38 + t36 * t52, -t53, 0, 0; 0, t50 + ((t34 * r_i_i_C(1) + t35 * r_i_i_C(2) + qJ(3)) * t38 + (-t35 * r_i_i_C(1) + t34 * r_i_i_C(2) - t55) * t36) * qJD(2), qJD(2) * t36, 0, 0, 0;];
JaD_transl  = t1;
