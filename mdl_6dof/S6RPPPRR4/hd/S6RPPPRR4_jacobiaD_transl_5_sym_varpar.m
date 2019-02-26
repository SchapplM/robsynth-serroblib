% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:17
% EndTime: 2019-02-26 20:24:18
% DurationCPUTime: 0.09s
% Computational Cost: add. (59->26), mult. (162->42), div. (0->0), fcn. (140->6), ass. (0->20)
t61 = -pkin(1) - pkin(2);
t60 = cos(pkin(9));
t50 = cos(qJ(1));
t59 = qJD(1) * t50;
t47 = sin(qJ(5));
t58 = qJD(5) * t47;
t49 = cos(qJ(5));
t57 = qJD(5) * t49;
t56 = pkin(3) + pkin(7) + r_i_i_C(3);
t48 = sin(qJ(1));
t55 = t48 * t60;
t54 = -r_i_i_C(1) * t47 - r_i_i_C(2) * t49 - qJ(4);
t53 = (r_i_i_C(1) * t49 - r_i_i_C(2) * t47) * qJD(5);
t46 = sin(pkin(9));
t52 = t48 * t46 + t50 * t60;
t51 = qJD(4) + t53;
t42 = t50 * t46 - t55;
t40 = -qJD(1) * t55 + t46 * t59;
t39 = t52 * qJD(1);
t1 = [t50 * qJD(2) - t51 * t52 + t54 * t40 - t56 * t39 + (-qJ(2) * t48 + t61 * t50) * qJD(1), t59, 0, -t39 (t39 * t47 - t42 * t57) * r_i_i_C(2) + (-t39 * t49 - t42 * t58) * r_i_i_C(1), 0; t48 * qJD(2) + t51 * t42 + t56 * t40 + t54 * t39 + (qJ(2) * t50 + t61 * t48) * qJD(1), qJD(1) * t48, 0, t40 (-t40 * t47 - t52 * t57) * r_i_i_C(2) + (t40 * t49 - t52 * t58) * r_i_i_C(1), 0; 0, 0, 0, 0, t53, 0;];
JaD_transl  = t1;
