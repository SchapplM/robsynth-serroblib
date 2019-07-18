% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR3_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_3_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:38
% DurationCPUTime: 0.10s
% Computational Cost: add. (75->24), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
t34 = qJD(2) + qJD(3);
t35 = qJ(2) + qJ(3);
t33 = cos(t35);
t53 = r_i_i_C(2) * t33;
t32 = sin(t35);
t55 = r_i_i_C(1) * t32;
t44 = t53 + t55;
t43 = t44 * t34;
t36 = sin(qJ(2));
t56 = pkin(1) * t36;
t57 = qJD(2) * t56 + t43;
t54 = r_i_i_C(2) * t32;
t52 = t33 * t34;
t37 = sin(qJ(1));
t51 = qJD(1) * t37;
t39 = cos(qJ(1));
t50 = qJD(1) * t39;
t38 = cos(qJ(2));
t49 = qJD(2) * t38;
t48 = r_i_i_C(1) * t52;
t47 = t34 * t54;
t46 = qJD(1) * t53;
t42 = -pkin(1) * t38 - r_i_i_C(1) * t33 + t54;
t41 = t37 * t46 + t51 * t55 + (t47 - t48) * t39;
t28 = t37 * t47;
t1 = [t57 * t37 + (-r_i_i_C(3) * t37 + t42 * t39) * qJD(1), (t36 * t51 - t39 * t49) * pkin(1) + t41, t41, 0, 0; -t57 * t39 + (r_i_i_C(3) * t39 + t42 * t37) * qJD(1), t28 + (-pkin(1) * t49 - t48) * t37 + (-t44 - t56) * t50, -t39 * t46 + t28 + (-t32 * t50 - t37 * t52) * r_i_i_C(1), 0, 0; 0, -t57, -t43, 0, 0;];
JaD_transl  = t1;
