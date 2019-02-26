% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:35
% EndTime: 2019-02-26 20:22:35
% DurationCPUTime: 0.08s
% Computational Cost: add. (64->23), mult. (84->35), div. (0->0), fcn. (54->6), ass. (0->16)
t25 = sin(qJ(5));
t26 = cos(qJ(5));
t28 = (r_i_i_C(1) * t26 - r_i_i_C(2) * t25) * qJD(5);
t36 = qJD(4) + t28;
t24 = qJ(1) + pkin(9);
t22 = sin(t24);
t35 = qJD(1) * t22;
t34 = qJD(1) * t25;
t33 = qJD(1) * t26;
t32 = qJD(5) * t25;
t31 = qJD(5) * t26;
t30 = pkin(7) + r_i_i_C(3) - qJ(3);
t27 = -r_i_i_C(1) * t25 - r_i_i_C(2) * t26 - pkin(2) - qJ(4);
t23 = cos(t24);
t21 = qJD(1) * t23;
t1 = [t23 * qJD(3) - t36 * t22 + (-cos(qJ(1)) * pkin(1) + t30 * t22 + t27 * t23) * qJD(1), 0, t21, -t35 (t22 * t34 - t23 * t31) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0; t22 * qJD(3) + t36 * t23 + (-sin(qJ(1)) * pkin(1) - t30 * t23 + t27 * t22) * qJD(1), 0, t35, t21 (-t22 * t31 - t23 * t34) * r_i_i_C(2) + (-t22 * t32 + t23 * t33) * r_i_i_C(1), 0; 0, 0, 0, 0, -t28, 0;];
JaD_transl  = t1;
