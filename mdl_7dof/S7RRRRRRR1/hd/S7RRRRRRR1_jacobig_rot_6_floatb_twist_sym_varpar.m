% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Jg_rot [3x7]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S7RRRRRRR1_jacobig_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobig_rot_6_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobig_rot_6_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:54
% EndTime: 2018-11-26 21:20:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->17), mult. (52->31), div. (0->0), fcn. (83->10), ass. (0->22)
t175 = sin(qJ(3));
t176 = sin(qJ(2));
t189 = t176 * t175;
t180 = cos(qJ(3));
t188 = t176 * t180;
t177 = sin(qJ(1));
t187 = t177 * t176;
t181 = cos(qJ(2));
t186 = t177 * t181;
t182 = cos(qJ(1));
t185 = t182 * t175;
t184 = t182 * t176;
t183 = t182 * t180;
t179 = cos(qJ(4));
t178 = cos(qJ(5));
t174 = sin(qJ(4));
t173 = sin(qJ(5));
t172 = -t177 * t175 + t181 * t183;
t171 = -t177 * t180 - t181 * t185;
t170 = t180 * t186 + t185;
t169 = -t175 * t186 + t183;
t1 = [0, t177, -t184, t171, t172 * t174 - t179 * t184 (t172 * t179 + t174 * t184) * t173 - t171 * t178, 0; 0, -t182, -t187, t169, t170 * t174 - t179 * t187 (t170 * t179 + t174 * t187) * t173 - t169 * t178, 0; 1, 0, t181, -t189, t174 * t188 + t181 * t179 (-t181 * t174 + t179 * t188) * t173 + t178 * t189, 0;];
Jg_rot  = t1;
