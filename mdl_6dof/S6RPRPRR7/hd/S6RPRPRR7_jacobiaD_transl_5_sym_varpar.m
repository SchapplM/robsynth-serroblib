% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:12
% EndTime: 2019-02-26 20:52:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (136->32), mult. (144->44), div. (0->0), fcn. (94->8), ass. (0->28)
t51 = qJ(3) + pkin(10);
t41 = cos(qJ(3)) * pkin(3) + pkin(4) * cos(t51);
t70 = t41 * qJD(3);
t47 = qJ(5) + t51;
t44 = cos(t47);
t69 = r_i_i_C(1) * t44;
t43 = sin(t47);
t68 = r_i_i_C(2) * t43;
t50 = qJD(3) + qJD(5);
t67 = t43 * t50;
t66 = t44 * t50;
t53 = sin(qJ(1));
t65 = qJD(1) * t53;
t55 = cos(qJ(1));
t48 = qJD(1) * t55;
t64 = -pkin(1) - r_i_i_C(3) - pkin(8) - qJ(4) - pkin(7);
t63 = r_i_i_C(1) * t67;
t61 = qJD(1) * t69;
t62 = t53 * t61 + (t66 * r_i_i_C(2) + t63) * t55;
t60 = -r_i_i_C(1) * t66 + r_i_i_C(2) * t67;
t40 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t51);
t59 = -r_i_i_C(1) * t43 - r_i_i_C(2) * t44;
t58 = qJD(1) * (t41 - t68);
t57 = qJ(2) + t40 - t59;
t56 = qJD(2) + t70 + (-t68 + t69) * t50;
t39 = t55 * t61;
t34 = t40 * qJD(3);
t1 = [-t53 * qJD(4) + t56 * t55 + (-t57 * t53 + t64 * t55) * qJD(1), t48, t39 + t55 * t58 + (t59 * t50 - t34) * t53, -t65, -t53 * t63 + t39 + (-t43 * t48 - t53 * t66) * r_i_i_C(2), 0; t55 * qJD(4) + t56 * t53 + (t64 * t53 + t57 * t55) * qJD(1), t65, t55 * t34 + t53 * t58 + t62, t48, -t65 * t68 + t62, 0; 0, 0, t60 - t70, 0, t60, 0;];
JaD_transl  = t1;
