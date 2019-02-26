% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JgD_rot = S6PRRRPP3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:02
% EndTime: 2019-02-26 20:10:02
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t151 = sin(pkin(6));
t154 = sin(qJ(3));
t164 = t151 * t154;
t155 = sin(qJ(2));
t163 = t151 * t155;
t153 = cos(pkin(6));
t162 = t153 * t155;
t157 = cos(qJ(2));
t161 = t153 * t157;
t160 = qJD(2) * t154;
t150 = sin(pkin(10));
t152 = cos(pkin(10));
t159 = t150 * t157 + t152 * t162;
t158 = -t150 * t162 + t152 * t157;
t156 = cos(qJ(3));
t1 = [0, 0, t158 * qJD(2) (t150 * t164 + t158 * t156) * qJD(3) + (-t150 * t161 - t152 * t155) * t160, 0, 0; 0, 0, t159 * qJD(2) (-t152 * t164 + t159 * t156) * qJD(3) + (-t150 * t155 + t152 * t161) * t160, 0, 0; 0, 0, qJD(2) * t163, t151 * t157 * t160 + (t153 * t154 + t156 * t163) * qJD(3), 0, 0;];
JgD_rot  = t1;
