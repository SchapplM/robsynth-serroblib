% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:04
% EndTime: 2019-02-26 21:13:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (56->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t83 = sin(qJ(3));
t84 = sin(qJ(1));
t91 = t84 * t83;
t82 = qJ(4) + qJ(5);
t80 = sin(t82);
t85 = cos(qJ(3));
t90 = t85 * t80;
t81 = cos(t82);
t89 = t85 * t81;
t86 = cos(qJ(1));
t88 = t86 * t83;
t87 = t86 * t85;
t79 = -t84 * t80 + t81 * t88;
t78 = t80 * t88 + t84 * t81;
t77 = t86 * t80 + t81 * t91;
t76 = -t80 * t91 + t86 * t81;
t1 = [t79, 0, t84 * t89, t76, t76, 0; t77, 0, -t81 * t87, t78, t78, 0; 0, 0, -t83 * t81, -t90, -t90, 0; -t78, 0, -t84 * t90, -t77, -t77, 0; t76, 0, t80 * t87, t79, t79, 0; 0, 0, t83 * t80, -t89, -t89, 0; -t87, 0, t91, 0, 0, 0; -t84 * t85, 0, -t88, 0, 0, 0; 0, 0, t85, 0, 0, 0;];
JR_rot  = t1;
