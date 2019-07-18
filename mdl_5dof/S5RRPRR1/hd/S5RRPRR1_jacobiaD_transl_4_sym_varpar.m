% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:52
% EndTime: 2019-07-18 17:22:52
% DurationCPUTime: 0.10s
% Computational Cost: add. (90->28), mult. (118->39), div. (0->0), fcn. (77->6), ass. (0->29)
t40 = sin(qJ(2));
t37 = qJD(2) + qJD(4);
t38 = qJ(2) + qJ(4);
t36 = cos(t38);
t61 = r_i_i_C(2) * t36;
t35 = sin(t38);
t63 = r_i_i_C(1) * t35;
t49 = t61 + t63;
t48 = t49 * t37;
t44 = pkin(2) + pkin(1);
t55 = qJD(2) * t44;
t64 = t40 * t55 + t48;
t62 = r_i_i_C(2) * t35;
t60 = r_i_i_C(3) + pkin(3) + qJ(3);
t59 = t36 * t37;
t58 = t40 * t44;
t41 = sin(qJ(1));
t57 = qJD(1) * t41;
t43 = cos(qJ(1));
t56 = qJD(1) * t43;
t54 = r_i_i_C(1) * t59;
t53 = t37 * t62;
t51 = qJD(1) * t61;
t52 = t41 * t51 + t43 * t53 + t57 * t63;
t42 = cos(qJ(2));
t47 = -t42 * t55 - t54;
t46 = -r_i_i_C(1) * t36 - t42 * t44 + t62;
t31 = t41 * t53;
t1 = [t43 * qJD(3) + t64 * t41 + (-t60 * t41 + t46 * t43) * qJD(1), t43 * t47 + t57 * t58 + t52, t56, -t43 * t54 + t52, 0; t41 * qJD(3) - t64 * t43 + (t46 * t41 + t60 * t43) * qJD(1), t31 + t47 * t41 + (-t49 - t58) * t56, t57, -t43 * t51 + t31 + (-t35 * t56 - t41 * t59) * r_i_i_C(1), 0; 0, -t64, 0, -t48, 0;];
JaD_transl  = t1;
