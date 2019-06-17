% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RPRRR1
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
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-12 14:37
% Revision: aab8d7cd0cba739f5e0ec8d53b8419901d1154b0 (2019-06-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-12 14:37:34
% EndTime: 2019-06-12 14:37:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (65->35), mult. (218->72), div. (0->0), fcn. (179->6), ass. (0->28)
t50 = cos(qJ(4));
t52 = cos(qJ(1));
t69 = t50 * t52;
t49 = sin(qJ(1));
t68 = qJD(1) * t49;
t67 = qJD(1) * t52;
t48 = sin(qJ(3));
t66 = qJD(3) * t48;
t51 = cos(qJ(3));
t65 = qJD(3) * t51;
t64 = qJD(3) * t52;
t63 = qJD(4) * t48;
t62 = qJD(4) * t51;
t61 = qJ(2) * qJD(1);
t60 = -qJD(1) + t62;
t59 = qJD(1) * t51 - qJD(4);
t47 = sin(qJ(4));
t58 = r_i_i_C(1) * t50 - r_i_i_C(2) * t47;
t57 = r_i_i_C(1) * t47 + r_i_i_C(2) * t50;
t56 = t60 * t47;
t55 = qJD(3) * t58;
t54 = -r_i_i_C(3) * qJD(3) + t57 * qJD(4);
t53 = t48 * t64 + t59 * t49;
t46 = -t59 * t69 + (t50 * t66 + t56) * t49;
t45 = t60 * t50 * t49 + (-t49 * t66 + t59 * t52) * t47;
t44 = t53 * t50 + t52 * t56;
t43 = t53 * t47 - t60 * t69;
t1 = [-t49 * t61 + t46 * r_i_i_C(1) + t45 * r_i_i_C(2) + t52 * qJD(2) + (-t48 * t67 - t49 * t65) * r_i_i_C(3), t67, (-r_i_i_C(3) * t68 - t52 * t55) * t51 + (t54 * t52 + t58 * t68) * t48, t43 * r_i_i_C(1) + t44 * r_i_i_C(2), 0; t52 * t61 - t44 * r_i_i_C(1) + t43 * r_i_i_C(2) + t49 * qJD(2) + (-t48 * t68 + t51 * t64) * r_i_i_C(3), t68, (r_i_i_C(3) * t67 - t49 * t55) * t51 + (t54 * t49 - t58 * t67) * t48, -t45 * r_i_i_C(1) + t46 * r_i_i_C(2), 0; 0, 0, -t57 * t62 + (r_i_i_C(3) * t51 - t58 * t48) * qJD(3), (t47 * t63 - t50 * t65) * r_i_i_C(2) + (-t47 * t65 - t50 * t63) * r_i_i_C(1), 0;];
JaD_transl  = t1;
