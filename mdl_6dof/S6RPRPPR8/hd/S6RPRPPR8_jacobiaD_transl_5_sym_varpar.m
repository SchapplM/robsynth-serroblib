% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiaD_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (62->24), mult. (176->36), div. (0->0), fcn. (121->4), ass. (0->17)
t20 = sin(qJ(3));
t22 = cos(qJ(3));
t27 = pkin(3) + pkin(4) - r_i_i_C(2);
t32 = r_i_i_C(1) + qJ(4);
t37 = -(t32 * t20 + t27 * t22) * qJD(3) + t22 * qJD(4);
t36 = -qJD(5) + (t27 * t20 - t32 * t22 + qJ(2)) * qJD(1);
t35 = t27 * qJD(3) - qJD(4);
t21 = sin(qJ(1));
t19 = qJD(1) * t21;
t23 = cos(qJ(1));
t31 = qJD(1) * t23;
t30 = qJD(3) * t21;
t29 = qJD(3) * t23;
t26 = qJD(1) * t32;
t25 = qJD(1) * t27;
t24 = qJD(2) + (-pkin(1) - pkin(7) + r_i_i_C(3) + qJ(5)) * qJD(1) - t37;
t1 = [-t36 * t21 + t24 * t23, t31 (t23 * t25 + t32 * t30) * t22 + (-t35 * t21 + t23 * t26) * t20, t20 * t30 - t22 * t31, t19, 0; t24 * t21 + t36 * t23, t19 (t21 * t25 - t32 * t29) * t22 + (t21 * t26 + t35 * t23) * t20, -t22 * t19 - t20 * t29, -t31, 0; 0, 0, t37, qJD(3) * t22, 0, 0;];
JaD_transl  = t1;
