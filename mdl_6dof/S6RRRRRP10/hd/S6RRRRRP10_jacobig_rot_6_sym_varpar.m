% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRP10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:01
% EndTime: 2019-02-26 22:45:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (14->9), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->18)
t158 = sin(pkin(6));
t162 = sin(qJ(1));
t171 = t162 * t158;
t161 = sin(qJ(2));
t170 = t162 * t161;
t164 = cos(qJ(2));
t169 = t162 * t164;
t165 = cos(qJ(1));
t168 = t165 * t158;
t167 = t165 * t161;
t166 = t165 * t164;
t163 = cos(qJ(3));
t160 = sin(qJ(3));
t159 = cos(pkin(6));
t157 = t158 * t161 * t160 - t159 * t163;
t156 = (-t159 * t170 + t166) * t160 - t163 * t171;
t155 = (t159 * t167 + t169) * t160 + t163 * t168;
t1 = [0, t171, t159 * t169 + t167, t156, t156, 0; 0, -t168, -t159 * t166 + t170, t155, t155, 0; 1, t159, -t158 * t164, t157, t157, 0;];
Jg_rot  = t1;
