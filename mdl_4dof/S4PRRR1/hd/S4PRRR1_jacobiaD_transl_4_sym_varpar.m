% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:03
% EndTime: 2019-02-26 19:29:03
% DurationCPUTime: 0.05s
% Computational Cost: add. (86->11), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->14)
t51 = qJD(2) + qJD(3);
t57 = pkin(3) * t51;
t56 = pkin(2) * qJD(2);
t50 = pkin(7) + qJ(2);
t49 = qJ(3) + t50;
t47 = qJ(4) + t49;
t43 = sin(t47);
t44 = cos(t47);
t48 = qJD(4) + t51;
t55 = (-r_i_i_C(1) * t44 + r_i_i_C(2) * t43) * t48;
t54 = (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * t48;
t53 = -cos(t49) * t57 + t55;
t52 = -sin(t49) * t57 + t54;
t1 = [0, -cos(t50) * t56 + t53, t53, t55; 0, -sin(t50) * t56 + t52, t52, t54; 0, 0, 0, 0;];
JaD_transl  = t1;
