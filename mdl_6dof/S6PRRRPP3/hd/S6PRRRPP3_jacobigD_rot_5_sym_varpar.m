% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPP3_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobigD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:02
% EndTime: 2019-02-26 20:10:02
% DurationCPUTime: 0.08s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t150 = sin(pkin(6));
t153 = sin(qJ(3));
t163 = t150 * t153;
t154 = sin(qJ(2));
t162 = t150 * t154;
t152 = cos(pkin(6));
t161 = t152 * t154;
t156 = cos(qJ(2));
t160 = t152 * t156;
t159 = qJD(2) * t153;
t149 = sin(pkin(10));
t151 = cos(pkin(10));
t158 = t149 * t156 + t151 * t161;
t157 = -t149 * t161 + t151 * t156;
t155 = cos(qJ(3));
t1 = [0, 0, t157 * qJD(2) (t149 * t163 + t155 * t157) * qJD(3) + (-t149 * t160 - t151 * t154) * t159, 0, 0; 0, 0, t158 * qJD(2) (-t151 * t163 + t155 * t158) * qJD(3) + (-t149 * t154 + t151 * t160) * t159, 0, 0; 0, 0, qJD(2) * t162, t150 * t156 * t159 + (t152 * t153 + t155 * t162) * qJD(3), 0, 0;];
JgD_rot  = t1;
