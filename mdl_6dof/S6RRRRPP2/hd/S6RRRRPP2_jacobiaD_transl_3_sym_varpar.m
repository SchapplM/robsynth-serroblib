% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP2_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP2_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiaD_transl_3_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:31
% EndTime: 2019-02-26 22:25:31
% DurationCPUTime: 0.09s
% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
t38 = qJD(2) + qJD(3);
t39 = qJ(2) + qJ(3);
t37 = cos(t39);
t59 = r_i_i_C(2) * t37;
t36 = sin(t39);
t61 = r_i_i_C(1) * t36;
t49 = t59 + t61;
t47 = t49 * t38;
t40 = sin(qJ(2));
t62 = pkin(2) * t40;
t63 = qJD(2) * t62 + t47;
t60 = r_i_i_C(2) * t36;
t58 = r_i_i_C(3) + pkin(8) + pkin(7);
t57 = t37 * t38;
t41 = sin(qJ(1));
t56 = qJD(1) * t41;
t43 = cos(qJ(1));
t55 = qJD(1) * t43;
t42 = cos(qJ(2));
t54 = qJD(2) * t42;
t53 = r_i_i_C(1) * t57;
t52 = t38 * t60;
t51 = qJD(1) * t59;
t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
t31 = t41 * t52;
t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1) (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
JaD_transl  = t1;
