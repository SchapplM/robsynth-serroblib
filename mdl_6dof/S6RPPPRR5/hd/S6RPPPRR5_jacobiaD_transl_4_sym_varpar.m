% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:53
% EndTime: 2019-02-26 20:24:53
% DurationCPUTime: 0.05s
% Computational Cost: add. (19->15), mult. (44->22), div. (0->0), fcn. (32->4), ass. (0->11)
t55 = -pkin(1) - qJ(3);
t54 = pkin(3) + qJ(2);
t51 = sin(qJ(1));
t53 = qJD(1) * t51;
t52 = cos(qJ(1));
t50 = cos(pkin(9));
t49 = sin(pkin(9));
t48 = qJD(1) * t52;
t47 = (-t49 * t51 + t50 * t52) * qJD(1);
t46 = (-t49 * t52 - t50 * t51) * qJD(1);
t1 = [t46 * r_i_i_C(1) - t47 * r_i_i_C(2) + t52 * qJD(2) - t51 * qJD(3) + (-t54 * t51 + t55 * t52) * qJD(1), t48, -t53, 0, 0, 0; t47 * r_i_i_C(1) + t46 * r_i_i_C(2) + t51 * qJD(2) + t52 * qJD(3) + (t55 * t51 + t54 * t52) * qJD(1), t53, t48, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JaD_transl  = t1;
