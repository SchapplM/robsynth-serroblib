% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR4_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR4_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:27:24
% EndTime: 2019-02-26 20:27:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (49->23), mult. (138->41), div. (0->0), fcn. (116->6), ass. (0->17)
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t51 = (r_i_i_C(1) * t47 + r_i_i_C(2) * t49) * qJD(4);
t57 = -pkin(1) - pkin(2);
t56 = -pkin(7) - r_i_i_C(3);
t55 = qJD(4) * t47;
t54 = qJD(4) * t49;
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t42 = t45 * t50 - t48 * t46;
t43 = t48 * t45 + t46 * t50;
t52 = r_i_i_C(1) * t49 - r_i_i_C(2) * t47 + pkin(3);
t41 = t42 * qJD(1);
t40 = t43 * qJD(1);
t1 = [t50 * qJD(2) + t56 * t41 - t42 * t51 - t52 * t40 + (-qJ(2) * t48 + t57 * t50) * qJD(1), qJD(1) * t50, 0 (-t41 * t49 + t43 * t55) * r_i_i_C(2) + (-t41 * t47 - t43 * t54) * r_i_i_C(1), 0, 0; t48 * qJD(2) + t56 * t40 - t43 * t51 + t52 * t41 + (qJ(2) * t50 + t57 * t48) * qJD(1), qJD(1) * t48, 0 (-t40 * t49 - t42 * t55) * r_i_i_C(2) + (-t40 * t47 + t42 * t54) * r_i_i_C(1), 0, 0; 0, 0, 0, t51, 0, 0;];
JaD_transl  = t1;
