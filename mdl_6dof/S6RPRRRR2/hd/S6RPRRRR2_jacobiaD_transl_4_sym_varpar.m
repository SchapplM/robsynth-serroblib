% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:35
% EndTime: 2019-02-26 21:15:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (119->29), mult. (118->40), div. (0->0), fcn. (75->8), ass. (0->26)
t44 = qJD(3) + qJD(4);
t46 = qJ(3) + qJ(4);
t42 = sin(t46);
t43 = cos(t46);
t63 = r_i_i_C(2) * t43;
t54 = r_i_i_C(1) * t42 + t63;
t52 = t54 * t44;
t47 = sin(qJ(3));
t65 = pkin(3) * t47;
t66 = qJD(3) * t65 + t52;
t64 = r_i_i_C(2) * t42;
t62 = r_i_i_C(3) + pkin(8) + pkin(7);
t61 = t43 * t44;
t60 = qJD(1) * t42;
t48 = cos(qJ(3));
t59 = qJD(3) * t48;
t58 = r_i_i_C(1) * t61;
t57 = t44 * t64;
t56 = qJD(1) * t63;
t53 = -pkin(3) * t48 - r_i_i_C(1) * t43 - pkin(2) + t64;
t45 = qJ(1) + pkin(11);
t40 = sin(t45);
t41 = cos(t45);
t51 = (t57 - t58) * t41 + (t60 * r_i_i_C(1) + t56) * t40;
t35 = t40 * t57;
t1 = [t66 * t40 + (-cos(qJ(1)) * pkin(1) - t62 * t40 + t53 * t41) * qJD(1), 0 (qJD(1) * t40 * t47 - t41 * t59) * pkin(3) + t51, t51, 0, 0; -t66 * t41 + (-sin(qJ(1)) * pkin(1) + t62 * t41 + t53 * t40) * qJD(1), 0, t35 + (-pkin(3) * t59 - t58) * t40 + (-t54 - t65) * t41 * qJD(1), -t41 * t56 + t35 + (-t40 * t61 - t41 * t60) * r_i_i_C(1), 0, 0; 0, 0, -t66, -t52, 0, 0;];
JaD_transl  = t1;
