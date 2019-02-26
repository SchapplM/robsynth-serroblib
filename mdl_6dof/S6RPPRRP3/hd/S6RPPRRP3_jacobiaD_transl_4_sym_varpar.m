% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:22
% EndTime: 2019-02-26 20:31:22
% DurationCPUTime: 0.10s
% Computational Cost: add. (53->19), mult. (76->33), div. (0->0), fcn. (48->6), ass. (0->14)
t24 = sin(qJ(4));
t25 = cos(qJ(4));
t34 = (r_i_i_C(1) * t25 - r_i_i_C(2) * t24) * qJD(4);
t33 = qJD(1) * t24;
t32 = qJD(1) * t25;
t31 = qJD(4) * t24;
t30 = qJD(4) * t25;
t29 = -pkin(2) - pkin(7) - r_i_i_C(3);
t27 = r_i_i_C(1) * t24 + r_i_i_C(2) * t25 + qJ(3);
t26 = qJD(3) + t34;
t23 = qJ(1) + pkin(9);
t22 = cos(t23);
t21 = sin(t23);
t1 = [t26 * t22 + (-cos(qJ(1)) * pkin(1) + t29 * t22 - t27 * t21) * qJD(1), 0, qJD(1) * t22 (-t21 * t30 - t22 * t33) * r_i_i_C(2) + (-t21 * t31 + t22 * t32) * r_i_i_C(1), 0, 0; t26 * t21 + (-sin(qJ(1)) * pkin(1) + t27 * t22 + t29 * t21) * qJD(1), 0, qJD(1) * t21 (-t21 * t33 + t22 * t30) * r_i_i_C(2) + (t21 * t32 + t22 * t31) * r_i_i_C(1), 0, 0; 0, 0, 0, -t34, 0, 0;];
JaD_transl  = t1;
