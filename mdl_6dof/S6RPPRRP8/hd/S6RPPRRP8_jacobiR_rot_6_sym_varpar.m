% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:12
% EndTime: 2019-02-26 20:34:12
% DurationCPUTime: 0.03s
% Computational Cost: add. (37->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t84 = sin(qJ(5));
t85 = sin(qJ(1));
t91 = t85 * t84;
t86 = cos(qJ(5));
t90 = t85 * t86;
t87 = cos(qJ(1));
t89 = t87 * t84;
t88 = t87 * t86;
t83 = pkin(9) + qJ(4);
t82 = cos(t83);
t81 = sin(t83);
t80 = t81 * t88 - t91;
t79 = t81 * t89 + t90;
t78 = t81 * t90 + t89;
t77 = t81 * t91 - t88;
t1 = [t80, 0, 0, t82 * t90, -t77, 0; t78, 0, 0, -t82 * t88, t79, 0; 0, 0, 0, -t81 * t86, -t82 * t84, 0; -t87 * t82, 0, 0, t85 * t81, 0, 0; -t85 * t82, 0, 0, -t87 * t81, 0, 0; 0, 0, 0, t82, 0, 0; t79, 0, 0, t82 * t91, t78, 0; t77, 0, 0, -t82 * t89, -t80, 0; 0, 0, 0, -t81 * t84, t82 * t86, 0;];
JR_rot  = t1;
