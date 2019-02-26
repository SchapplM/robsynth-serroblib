% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:47
% EndTime: 2019-02-26 20:26:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (97->21), mult. (106->27), div. (0->0), fcn. (69->8), ass. (0->17)
t31 = qJ(4) + pkin(10);
t27 = sin(t31);
t29 = cos(t31);
t40 = sin(qJ(4)) * pkin(4) + r_i_i_C(1) * t27 + r_i_i_C(2) * t29;
t44 = qJD(4) * t40;
t39 = cos(qJ(4)) * pkin(4) + r_i_i_C(1) * t29 - r_i_i_C(2) * t27;
t43 = t39 * qJD(4);
t32 = qJ(1) + pkin(9);
t28 = sin(t32);
t42 = qJD(1) * t28;
t41 = -pkin(2) - r_i_i_C(3) - qJ(5) - pkin(7);
t38 = qJ(3) + t40;
t37 = qJD(1) * t39;
t36 = qJD(3) + t43;
t30 = cos(t32);
t26 = qJD(1) * t30;
t1 = [-t28 * qJD(5) + t36 * t30 + (-cos(qJ(1)) * pkin(1) + t41 * t30 - t38 * t28) * qJD(1), 0, t26, -t28 * t44 + t30 * t37, -t42, 0; t30 * qJD(5) + t36 * t28 + (-sin(qJ(1)) * pkin(1) + t41 * t28 + t38 * t30) * qJD(1), 0, t42, t28 * t37 + t30 * t44, t26, 0; 0, 0, 0, -t43, 0, 0;];
JaD_transl  = t1;
