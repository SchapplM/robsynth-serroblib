% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR4_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR4_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:27
% EndTime: 2019-02-26 20:24:27
% DurationCPUTime: 0.06s
% Computational Cost: add. (28->17), mult. (70->23), div. (0->0), fcn. (58->4), ass. (0->13)
t31 = -pkin(1) - pkin(2);
t30 = pkin(3) - r_i_i_C(2);
t29 = -r_i_i_C(3) - qJ(4);
t24 = sin(qJ(1));
t28 = qJD(1) * t24;
t25 = cos(qJ(1));
t27 = qJD(1) * t25;
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t26 = t24 * t22 + t25 * t23;
t20 = t22 * t27 - t23 * t28;
t19 = t26 * qJD(1);
t1 = [-t26 * qJD(4) + t25 * qJD(2) + t29 * t20 - t30 * t19 + (-t24 * qJ(2) + t31 * t25) * qJD(1), t27, 0, -t19, 0, 0; -(-t25 * t22 + t24 * t23) * qJD(4) + t24 * qJD(2) + t30 * t20 + t29 * t19 + (t25 * qJ(2) + t31 * t24) * qJD(1), t28, 0, t20, 0, 0; 0, 0, 0, 0, 0, 0;];
JaD_transl  = t1;
