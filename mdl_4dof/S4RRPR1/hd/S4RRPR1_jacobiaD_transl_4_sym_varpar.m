% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRPR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_jacobiaD_transl_4_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:36:10
% EndTime: 2019-02-26 19:36:10
% DurationCPUTime: 0.06s
% Computational Cost: add. (96->13), mult. (44->15), div. (0->0), fcn. (22->8), ass. (0->13)
t60 = pkin(1) * qJD(1);
t55 = qJ(1) + qJ(2);
t54 = qJD(1) + qJD(2);
t51 = pkin(7) + t55;
t49 = qJ(4) + t51;
t45 = sin(t49);
t46 = cos(t49);
t50 = qJD(4) + t54;
t59 = (-r_i_i_C(1) * t46 + r_i_i_C(2) * t45) * t50;
t58 = (-r_i_i_C(1) * t45 - r_i_i_C(2) * t46) * t50;
t57 = (-pkin(2) * cos(t55) - pkin(3) * cos(t51)) * t54 + t59;
t56 = (-pkin(2) * sin(t55) - pkin(3) * sin(t51)) * t54 + t58;
t1 = [-cos(qJ(1)) * t60 + t57, t57, 0, t59; -sin(qJ(1)) * t60 + t56, t56, 0, t58; 0, 0, 0, 0;];
JaD_transl  = t1;
