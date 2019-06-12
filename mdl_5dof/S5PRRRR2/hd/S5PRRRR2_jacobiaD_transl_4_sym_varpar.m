% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5PRRRR2
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-03 15:11
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5PRRRR2_jacobiaD_transl_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-03 15:11:34
% EndTime: 2019-06-03 15:11:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (75->24), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
t36 = qJD(3) + qJD(4);
t37 = qJ(3) + qJ(4);
t35 = cos(t37);
t55 = r_i_i_C(2) * t35;
t34 = sin(t37);
t56 = r_i_i_C(1) * t34;
t44 = t55 + t56;
t38 = sin(qJ(3));
t57 = pkin(1) * t38;
t45 = qJD(3) * t57;
t58 = t36 * t44 + t45;
t54 = t34 * t36;
t53 = t35 * t36;
t52 = r_i_i_C(1) * t54 + r_i_i_C(2) * t53;
t39 = sin(qJ(2));
t51 = qJD(2) * t39;
t41 = cos(qJ(2));
t50 = qJD(2) * t41;
t40 = cos(qJ(3));
t49 = qJD(3) * t40;
t48 = r_i_i_C(1) * t53;
t47 = r_i_i_C(2) * t54;
t46 = qJD(2) * t55;
t43 = -pkin(1) * t40 - r_i_i_C(1) * t35 + r_i_i_C(2) * t34;
t42 = t39 * t46 + t51 * t56 + (t47 - t48) * t41;
t28 = t39 * t47;
t1 = [0, t58 * t39 + (-r_i_i_C(3) * t39 + t41 * t43) * qJD(2), (t38 * t51 - t41 * t49) * pkin(1) + t42, t42, 0; 0, 0, t45 + t52, t52, 0; 0, -t58 * t41 + (r_i_i_C(3) * t41 + t39 * t43) * qJD(2), t28 + (-pkin(1) * t49 - t48) * t39 + (-t44 - t57) * t50, -t41 * t46 + t28 + (-t34 * t50 - t39 * t53) * r_i_i_C(1), 0;];
JaD_transl  = t1;
