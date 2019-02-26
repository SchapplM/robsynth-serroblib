% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:47
% EndTime: 2019-02-26 20:23:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (32->17), mult. (84->24), div. (0->0), fcn. (70->6), ass. (0->12)
t36 = -pkin(1) - pkin(2);
t35 = -r_i_i_C(3) - qJ(4);
t27 = sin(pkin(9));
t29 = cos(pkin(9));
t30 = sin(qJ(1));
t31 = cos(qJ(1));
t34 = t31 * t27 - t30 * t29;
t33 = t30 * t27 + t31 * t29;
t32 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(3);
t25 = t34 * qJD(1);
t24 = t33 * qJD(1);
t1 = [-t33 * qJD(4) + t31 * qJD(2) + t35 * t25 - t32 * t24 + (-t30 * qJ(2) + t36 * t31) * qJD(1), qJD(1) * t31, 0, -t24, 0, 0; t34 * qJD(4) + t30 * qJD(2) + t35 * t24 + t32 * t25 + (t31 * qJ(2) + t36 * t30) * qJD(1), qJD(1) * t30, 0, t25, 0, 0; 0, 0, 0, 0, 0, 0;];
JaD_transl  = t1;
