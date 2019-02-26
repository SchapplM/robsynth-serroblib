% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:19
% EndTime: 2019-02-26 20:40:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (121->24), mult. (172->33), div. (0->0), fcn. (117->6), ass. (0->18)
t26 = sin(qJ(3));
t27 = cos(qJ(3));
t35 = pkin(3) + pkin(4) - r_i_i_C(2);
t41 = r_i_i_C(1) + qJ(4);
t32 = t35 * t26 - t41 * t27;
t42 = t32 * qJD(3) - t26 * qJD(4);
t25 = qJ(1) + pkin(9);
t23 = sin(t25);
t40 = qJD(1) * t23;
t24 = cos(t25);
t39 = qJD(1) * t24;
t38 = qJD(1) * t26;
t37 = qJD(3) * t27;
t34 = pkin(7) - r_i_i_C(3) - qJ(5);
t33 = -t41 * t26 - t35 * t27;
t30 = -pkin(2) + t33;
t29 = t33 * qJD(3) + qJD(4) * t27;
t1 = [-t24 * qJD(5) + t42 * t23 + (-cos(qJ(1)) * pkin(1) - t34 * t23 + t30 * t24) * qJD(1), 0, t29 * t24 + t32 * t40, -t23 * t38 + t24 * t37, -t39, 0; -t23 * qJD(5) - t42 * t24 + (-sin(qJ(1)) * pkin(1) + t34 * t24 + t30 * t23) * qJD(1), 0, t29 * t23 - t32 * t39, t23 * t37 + t24 * t38, -t40, 0; 0, 0, -t42, qJD(3) * t26, 0, 0;];
JaD_transl  = t1;
