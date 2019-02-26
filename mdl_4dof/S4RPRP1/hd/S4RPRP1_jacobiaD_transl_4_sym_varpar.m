% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPRP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPRP1_jacobiaD_transl_4_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_jacobiaD_transl_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_jacobiaD_transl_4_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPRP1_jacobiaD_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_jacobiaD_transl_4_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:45
% EndTime: 2018-11-14 13:48:45
% DurationCPUTime: 0.06s
% Computational Cost: add. (84->13), mult. (46->15), div. (0->0), fcn. (26->6), ass. (0->12)
t32 = qJ(4) + r_i_i_C(3);
t31 = -pkin(3) - r_i_i_C(1);
t26 = qJ(1) + pkin(6);
t24 = qJ(3) + t26;
t22 = sin(t24);
t25 = qJD(1) + qJD(3);
t30 = t25 * t22;
t23 = cos(t24);
t29 = t25 * t23;
t28 = t22 * qJD(4) + t29 * t32 + t31 * t30;
t27 = t23 * qJD(4) + (-t22 * t32 + t31 * t23) * t25;
t1 = [(-cos(t26) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t27, 0, t27, t29; (-sin(t26) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t28, 0, t28, t30; 0, 0, 0, 0;];
JaD_transl  = t1;
